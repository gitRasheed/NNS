// partial_moments.cpp
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>
#include "partial_moments.h"

using namespace Rcpp;
using namespace RcppParallel;

// Dispatching to the thread pool costs tens of microseconds per call, more than
// the arithmetic on small inputs (an 80 x 3 matrix takes about 25 us on the
// calling thread and 110 us dispatched). Below NNS_SERIAL_WORK element
// evaluations the worker runs on the calling thread instead. Every worker here
// writes only its own output slots, so the results are identical either way.
static const double NNS_SERIAL_WORK = 4096.0;

template <typename W>
static inline void nns_parallel_for(std::size_t begin, std::size_t end, W& worker, double work) {
  if (end <= begin) return;
  if (work < NNS_SERIAL_WORK) worker(begin, end);
  else parallelFor(begin, end, worker);
}

// A partial-moment total below this (or non-finite) has lost relative precision
// to underflow, or overflowed: terms below DBL_MIN are denormal or zero, so the
// total is only accurate to about n * DBL_MIN in absolute terms.
static const double NNS_SAFE_TOTAL = DBL_MIN / DBL_EPSILON;

static inline bool nns_unsafe_total(double total) {
  return !(total >= NNS_SAFE_TOTAL && total <= DBL_MAX);
}

static inline double lower_component(double diff, double degree, bool degree_is_int) {
  if (degree == 0) return diff >= 0.0 ? 1.0 : 0.0;
  if (diff < 0.0) return 0.0;
  return degree_is_int
  ? nns_int_pow(diff, static_cast<int>(degree))
    : std::pow(diff, degree);
}

static inline double upper_component(double diff, double degree, bool degree_is_int) {
  if (degree == 0) return diff > 0.0 ? 1.0 : 0.0;
  if (diff < 0.0) return 0.0;
  return degree_is_int
  ? nns_int_pow(diff, static_cast<int>(degree))
    : std::pow(diff, degree);
}

/////////////////
// UPM / LPM
// single thread
double LPM_C(const double &degree, const double &target, const RVector<double> &variable) {
  size_t n = variable.size();
  double out = 0;
  double value;
  
  for (size_t i = 0; i < n; i++) {
    value = target - variable[i];
    if (value >= 0) {
      if (nns_is_integer_degree(degree)) {
        if (degree == 0) {
          out += 1;
        } else if (degree == 1) {
          out += value;
        } else {
          out += nns_int_pow(value, static_cast<int>(degree));
        }
      } else {
        out += std::pow(value, degree);
      }
    } else out+= 0;
  }
  out /= n;
  return out;
}

double UPM_C(const double &degree, const double &target, const RVector<double> &variable) {
  size_t n = variable.size();
  double out = 0;
  double value;
  
  for (size_t i = 0; i < n; i++) {
    value = variable[i] - target;
    if (value > 0) {
      if (nns_is_integer_degree(degree)) {
        if (degree == 0) {
          out += 1;
        } else if (degree == 1) {
          out += value;
        } else {
          out += nns_int_pow(value, static_cast<int>(degree));
        }
      } else {
        out += std::pow(value, degree);
      }
    } else out+= 0;
  }
  out /= n;
  return out;
}


// Lower Partial Moment (LPM) count: degree == 0
struct CoLPM_CountWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  RVector<double> output;
  CoLPM_CountWorker(const NumericMatrix& data_, const NumericVector& target_, NumericVector& output_)
    : data(data_), target(target_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      bool below_all = true;
      for (std::size_t j = 0; j < d; ++j) {
        if (data(i, j) > target[j]) { below_all = false; break; }
      }
      output[i] = below_all ? 1.0 : 0.0;
    }
  }
};

// Lower Partial Moment (LPM) sum: degree > 0
struct CoLPM_SumWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  const double degree;
  RVector<double> output;
  CoLPM_SumWorker(const NumericMatrix& data_, const NumericVector& target_, double degree_, NumericVector& output_)
    : data(data_), target(target_), degree(degree_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      double prod = 1.0;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = target[j] - data(i, j);
        if (diff < 0.0) { prod = 0.0; break; }
        prod *= nns_is_integer_degree(degree)
          ? nns_int_pow(diff, static_cast<int>(degree))
            : std::pow(diff, degree);
      }
      output[i] = prod;
    }
  }
};

// Upper Partial Moment (UPM) count: degree == 0
struct CoUPM_CountWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  RVector<double> output;
  CoUPM_CountWorker(const NumericMatrix& data_, const NumericVector& target_, NumericVector& output_)
    : data(data_), target(target_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      // Strictly above, as Co.UPM does: ties go to the lower side so that at
      // degree 0 CLPM + CUPM + DPM partition the sample.
      bool above_all = true;
      for (std::size_t j = 0; j < d; ++j) {
        if (data(i, j) <= target[j]) { above_all = false; break; }
      }
      output[i] = above_all ? 1.0 : 0.0;
    }
  }
};

// Upper Partial Moment (UPM) sum: degree > 0
struct CoUPM_SumWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  const double degree;
  RVector<double> output;
  CoUPM_SumWorker(const NumericMatrix& data_, const NumericVector& target_, double degree_, NumericVector& output_)
    : data(data_), target(target_), degree(degree_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      double prod = 1.0;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = data(i, j) - target[j];
        if (diff < 0.0) { prod = 0.0; break; }
        prod *= nns_is_integer_degree(degree)
          ? nns_int_pow(diff, static_cast<int>(degree))
            : std::pow(diff, degree);
      }
      output[i] = prod;
    }
  }
};

// Discordant Partial Moment (DPM) count: degree == 0
struct DpmCountWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  RVector<double> output;
  DpmCountWorker(const NumericMatrix& data_, const NumericVector& target_, NumericVector& output_)
    : data(data_), target(target_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      // Discordant = neither all-lower (x <= target) nor all-upper (x > target).
      bool allBelow = true, allAbove = true;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = data(i, j) - target[j];
        if (diff > 0.0) allBelow = false;
        if (diff <= 0.0) allAbove = false;
        if (!allBelow && !allAbove) break;
      }
      output[i] = (!allBelow && !allAbove) ? 1.0 : 0.0;
    }
  }
};

// Discordant Partial Moment (DPM) sum: degree > 0
struct DpmSumWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  const double degree;
  RVector<double> output;
  DpmSumWorker(const NumericMatrix& data_, const NumericVector& target_, double degree_, NumericVector& output_)
    : data(data_), target(target_), degree(degree_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      bool allBelow = true, allAbove = true;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = data(i, j) - target[j];
        if (diff >= 0.0) allBelow = false;
        if (diff <= 0.0) allAbove = false;
        if (!allBelow && !allAbove) break;
      }
      if (allBelow || allAbove) { output[i] = 0.0; continue; }
      double prod = 1.0;
      for (std::size_t j = 0; j < d; ++j) {
        double abs_dev = std::abs(data(i, j) - target[j]);
        prod *= nns_is_integer_degree(degree)
          ? nns_int_pow(abs_dev, static_cast<int>(degree))
            : std::pow(abs_dev, degree);
      }
      output[i] = prod;
    }
  }
};

// Normalized nD partial moments (CLPM, CUPM, DPM) / (CLPM + CUPM + DPM) for a
// target whose plain total is unsafe (zero, under- or overflowed). Pure C++ with
// no R API, so the batch worker can call it.
//
// For degree > 0 a row contributes to none of the three cells iff it ties the
// target in some column. If every row ties, the total is exactly 0 (0/0): return
// the degree-0 cells (inclusive all-below, strict all-above, the rest), as
// LPM.ratio does for a constant variable; they partition the sample. Otherwise
// the zero or non-finite total is under/overflow: divide column j's deviations by
// s_j = max_i |x_ij - t_j| (the cells are invariant to per-column scale) and, if
// the products still underflow, sum them in the log domain.
static void nD_normalized_rescue(const RMatrix<double> &data, const double *t,
                                 const double degree, double out[3]) {
  const std::size_t n = data.nrow(), d = data.ncol();
  std::vector<double> s(d, 0.0);
  for (std::size_t j = 0; j < d; ++j)
    for (std::size_t i = 0; i < n; ++i)
      s[j] = std::max(s[j], std::fabs(data(i, j) - t[j]));
  
  // cell: 0 = CLPM (all below), 1 = CUPM (all above), 2 = DPM, -1 = tied row.
  std::vector<int> cell(n);
  std::vector<double> logterm(n);
  double direct[3] = {0.0, 0.0, 0.0};
  double c0[3] = {0.0, 0.0, 0.0};
  double max_log = -std::numeric_limits<double>::infinity();
  bool any_untied = false;
  for (std::size_t i = 0; i < n; ++i) {
    bool all_le = true, all_gt = true, all_lt = true, tied = false;
    double prod = 1.0, lsum = 0.0;
    for (std::size_t j = 0; j < d; ++j) {
      const double diff = data(i, j) - t[j];
      if (diff > 0.0) all_le = false;
      if (diff <= 0.0) all_gt = false;
      if (diff >= 0.0) all_lt = false;
      if (diff == 0.0) { tied = true; continue; }
      const double z = std::fabs(diff) / s[j];
      prod *= std::pow(z, degree);
      lsum += std::log(z);
    }
    c0[all_le ? 0 : (all_gt ? 1 : 2)] += 1.0;
    if (tied) { cell[i] = -1; continue; }
    any_untied = true;
    cell[i] = all_lt ? 0 : (all_gt ? 1 : 2);
    direct[cell[i]] += prod;
    logterm[i] = degree * lsum;
    max_log = std::max(max_log, logterm[i]);
  }
  
  if (!any_untied) {
    for (int k = 0; k < 3; ++k) out[k] = c0[k] / double(n);
    return;
  }
  double total = direct[0] + direct[1] + direct[2];
  if (nns_unsafe_total(total)) {
    direct[0] = direct[1] = direct[2] = 0.0;
    for (std::size_t i = 0; i < n; ++i)
      if (cell[i] >= 0) direct[cell[i]] += std::exp(logterm[i] - max_log);
    total = direct[0] + direct[1] + direct[2];
  }
  for (int k = 0; k < 3; ++k) out[k] = direct[k] / total;
}

static inline double nD_normalized_rescue_cell(const NumericMatrix &data,
                                               const NumericVector &target,
                                               const double degree, const int k) {
  double out[3];
  nD_normalized_rescue(RMatrix<double>(data), &target[0], degree, out);
  return out[k];
}

double clpm_nD_cpp(const NumericMatrix& data,
                   const NumericVector& target,
                   double degree,
                   bool norm) {
  size_t n = data.nrow();
  size_t d = data.ncol();
  if (static_cast<size_t>(target.size()) != d)
    stop("`target` length must match number of columns in `data`");
  
  if (degree == 0.0) {
    NumericVector counts(n);
    CoLPM_CountWorker countWorker(data, target, counts);
    nns_parallel_for(0, n, countWorker, double(n) * double(d));
    return sum(counts) / double(n);
  }
  
  NumericVector vals(n);
  CoLPM_SumWorker sumWorker(data, target, degree, vals);
  nns_parallel_for(0, n, sumWorker, double(n) * double(d));
  double clpm_un = sum(vals) / double(n);
  double result = clpm_un;
  
  if (norm) {
    double cupm_un = cupm_nD_cpp(data, target, degree, false);
    double dpm_un  = dpm_nD_cpp(data, target, degree, false);
    double norm_const = clpm_un + cupm_un + dpm_un;
    result = nns_unsafe_total(norm_const)
      ? nD_normalized_rescue_cell(data, target, degree, 0)
      : clpm_un / norm_const;
  }
  return result;
}

double cupm_nD_cpp(const NumericMatrix& data,
                   const NumericVector& target,
                   double degree,
                   bool norm) {
  size_t n = data.nrow();
  size_t d = data.ncol();
  if (static_cast<size_t>(target.size()) != d)
    stop("`target` length must match number of columns in `data`");
  
  if (degree == 0.0) {
    NumericVector counts(n);
    CoUPM_CountWorker countWorker(data, target, counts);
    nns_parallel_for(0, n, countWorker, double(n) * double(d));
    return sum(counts) / double(n);
  }
  
  NumericVector vals(n);
  CoUPM_SumWorker sumWorker(data, target, degree, vals);
  nns_parallel_for(0, n, sumWorker, double(n) * double(d));
  double cupm_un = sum(vals) / double(n);
  double result = cupm_un;
  
  if (norm) {
    double clpm_un = clpm_nD_cpp(data, target, degree, false);
    double dpm_un  = dpm_nD_cpp(data, target, degree, false);
    double norm_const = clpm_un + cupm_un + dpm_un;
    result = nns_unsafe_total(norm_const)
      ? nD_normalized_rescue_cell(data, target, degree, 1)
      : cupm_un / norm_const;
  }
  return result;
}

double dpm_nD_cpp(const NumericMatrix& data,
                  const NumericVector& target,
                  double degree,
                  bool norm) {
  size_t n = data.nrow();
  size_t d = data.ncol();
  if (static_cast<size_t>(target.size()) != d)
    stop("`target` length must match number of columns in `data`");
  
  if (degree == 0.0) {
    NumericVector counts(n);
    DpmCountWorker countWorker(data, target, counts);
    nns_parallel_for(0, n, countWorker, double(n) * double(d));
    return sum(counts) / double(n);
  }
  
  NumericVector vals(n);
  DpmSumWorker sumWorker(data, target, degree, vals);
  nns_parallel_for(0, n, sumWorker, double(n) * double(d));
  double dpm_un = sum(vals) / double(n);
  double result = dpm_un;
  
  if (norm) {
    double clpm_un = clpm_nD_cpp(data, target, degree, false);
    double cupm_un = cupm_nD_cpp(data, target, degree, false);
    double norm_const = clpm_un + cupm_un + dpm_un;
    result = nns_unsafe_total(norm_const)
      ? nD_normalized_rescue_cell(data, target, degree, 2)
      : dpm_un / norm_const;
  }
  return result;
}


// ============================================================================
// Batched nD CoLPM backend
// ============================================================================
//
// Computes CoLPM_nD(data, target_row, degree, norm) for every row of `targets`
// in one C++ call.  This replaces the R-side pattern:
//   apply(variable, 1, function(row) Co.LPM_nD(variable, row, degree = degree))
//
// Semantics match clpm_nD_cpp() (cupm_nD_cpp() with upper = true):
//   degree == 0 returns the raw lower-count probability, regardless of norm.
//   degree  > 0 returns raw CLPM if norm = false.
//   degree  > 0 returns CLPM / (CLPM + CUPM + DPM) if norm = true.

struct CoLPMnDBatchWorker : public Worker {
  const RMatrix<double> data;
  const RMatrix<double> targets;
  const double degree;
  const bool norm;
  const bool upper;
  const bool degree_is_int;
  RVector<double> output;
  
  CoLPMnDBatchWorker(const NumericMatrix& data_,
                     const NumericMatrix& targets_,
                     double degree_,
                     bool norm_,
                     bool upper_,
                     NumericVector& output_)
    : data(data_),
      targets(targets_),
      degree(degree_),
      norm(norm_),
      upper(upper_),
      degree_is_int(nns_is_integer_degree(degree_)),
      output(output_) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    const std::size_t n_obs = data.nrow();
    const std::size_t d = data.ncol();
    
    for (std::size_t r = begin; r < end; ++r) {
      
      // Match clpm_nD_cpp degree == 0 behavior:
      // it returns raw count probability and does not apply norm.
      if (degree == 0.0) {
        double count = 0.0;
        
        for (std::size_t i = 0; i < n_obs; ++i) {
          bool below_all = true;
          
          for (std::size_t j = 0; j < d; ++j) {
            // Inclusive below for CLPM, strictly above for CUPM.
            if (upper ? data(i, j) <= targets(r, j) : data(i, j) > targets(r, j)) {
              below_all = false;
              break;
            }
          }
          
          if (below_all) count += 1.0;
        }
        
        output[r] = count / static_cast<double>(n_obs);
        continue;
      }
      
      double clpm_sum = 0.0;
      double cupm_sum = 0.0;
      double dpm_sum  = 0.0;
      
      for (std::size_t i = 0; i < n_obs; ++i) {
        double lower_prod = 1.0;
        double upper_prod = 1.0;
        double dpm_prod   = 1.0;
        
        bool all_below_strict = true;
        bool all_above_strict = true;
        
        for (std::size_t j = 0; j < d; ++j) {
          const double diff = data(i, j) - targets(r, j);
          
          // CLPM component: target - data
          lower_prod *= lower_component(-diff, degree, degree_is_int);
          
          // CUPM component: data - target
          upper_prod *= upper_component(diff, degree, degree_is_int);
          
          // Match DpmSumWorker strict all-below/all-above logic.
          if (diff >= 0.0) all_below_strict = false;
          if (diff <= 0.0) all_above_strict = false;
          
          dpm_prod *= degree_is_int
          ? nns_int_pow(std::abs(diff), static_cast<int>(degree))
            : std::pow(std::abs(diff), degree);
        }
        
        clpm_sum += lower_prod;
        cupm_sum += upper_prod;
        
        if (!(all_below_strict || all_above_strict)) {
          dpm_sum += dpm_prod;
        }
      }
      
      const double inv_n = 1.0 / static_cast<double>(n_obs);
      const double clpm_un = clpm_sum * inv_n;
      const double cupm_un = cupm_sum * inv_n;
      const double num_un = upper ? cupm_un : clpm_un;
      
      if (!norm) {
        output[r] = num_un;
      } else {
        const double dpm_un  = dpm_sum  * inv_n;
        const double norm_const = clpm_un + cupm_un + dpm_un;
        
        if (nns_unsafe_total(norm_const)) {
          std::vector<double> t_row(d);
          for (std::size_t j = 0; j < d; ++j) t_row[j] = targets(r, j);
          double cells[3];
          nD_normalized_rescue(data, t_row.data(), degree, cells);
          output[r] = cells[upper ? 1 : 0];
        } else {
          output[r] = num_un / norm_const;
        }
      }
    }
  }
};


// [[Rcpp::export]]
NumericVector CoLPM_nD_batch_RCPP(const NumericMatrix& data,
                                  const NumericMatrix& targets,
                                  double degree = 0.0,
                                  bool norm = true,
                                  bool upper = false) {
  nns_check_degree(degree, "degree");
  if (data.ncol() != targets.ncol()) {
    stop("`targets` must have the same number of columns as `data`");
  }
  
  if (data.nrow() == 0) {
    stop("`data` must have at least one row");
  }
  
  NumericVector output(targets.nrow());
  
  CoLPMnDBatchWorker worker(data, targets, degree, norm, upper, output);
  nns_parallel_for(0, targets.nrow(), worker,
                   double(targets.nrow()) * double(data.nrow()) * double(data.ncol()));
  
  return output;
}

// parallelFor
#define NNS_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS)      \
size_t target_size=target.size();                        \
NumericVector output = NumericVector(target_size);       \
WORKER_CLASS tmp_func(degree, target, variable, output); \
nns_parallel_for(0, target_size, tmp_func, double(target_size) * double(variable.size())); \
return(output);

// Scalar guard: the prefix backend costs O(n log n + n*degree) to build,
// which only amortizes over many targets (crossover ~ log2(n) + degree).
// For few targets, a direct O(n) scan per target via the legacy kernels is
// faster and bit-identical to pre-13.0 semantics.
static const R_xlen_t NNS_DIRECT_PATH_MAX_TARGETS = 32;

// [[Rcpp::export]]
NumericVector LPM_CPv(const double &degree,
                      const NumericVector &target,
                      const NumericVector &variable) {
  nns_check_degree(degree, "degree");
  if (target.size() <= NNS_DIRECT_PATH_MAX_TARGETS) {
    NumericVector output(target.size());
    RcppParallel::RVector<double> v(variable);
    
    for (R_xlen_t i = 0; i < target.size(); ++i) {
      output[i] = LPM_C(degree, target[i], v);
    }
    
    return output;
  }
  
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(LPM_Worker);
}

// [[Rcpp::export]]
NumericVector UPM_CPv(const double &degree,
                      const NumericVector &target,
                      const NumericVector &variable) {
  nns_check_degree(degree, "degree");
  if (target.size() <= NNS_DIRECT_PATH_MAX_TARGETS) {
    NumericVector output(target.size());
    RcppParallel::RVector<double> v(variable);
    
    for (R_xlen_t i = 0; i < target.size(); ++i) {
      output[i] = UPM_C(degree, target[i], v);
    }
    
    return output;
  }
  
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(UPM_Worker);
}

// L/(L+U) is scale-free, but the raw sums of |x - t|^degree under- or overflow
// for small or large deviations at high degree (0/0 or Inf/Inf = NaN). For each
// such target, recompute with deviations divided by s = max|x - t|: the largest
// term is then exactly 1, so the total is at least 1 and finite. Targets with a
// safe total are left untouched; s == 0 (every observation equals the target)
// is left as computed, for fill_degenerate_ratio.
static NumericVector rescale_unsafe_ratio(NumericVector lpm, NumericVector upm,
                                          const double degree,
                                          const NumericVector &target,
                                          const NumericVector &variable,
                                          const bool lower) {
  const R_xlen_t n = variable.size();
  NumericVector out = (lower ? lpm : upm) / (lpm + upm);
  if (n == 0) return out;
  bool have_range = false;
  double vmin = 0.0, vmax = 0.0;
  for (R_xlen_t i = 0; i < out.size(); ++i) {
    if (!nns_unsafe_total(lpm[i] + upm[i])) continue;
    if (!have_range) {
      vmin = vmax = variable[0];
      for (R_xlen_t j = 1; j < n; ++j) {
        vmin = std::min(vmin, variable[j]);
        vmax = std::max(vmax, variable[j]);
      }
      have_range = true;
    }
    const double t = target[i];
    const double s = std::max(std::fabs(vmax - t), std::fabs(t - vmin));
    if (!(s > 0.0) || !R_FINITE(s)) continue;
    double L = 0.0, U = 0.0;
    for (R_xlen_t j = 0; j < n; ++j) {
      const double z = (variable[j] - t) / s;
      if (z < 0.0) L += std::pow(-z, degree);
      else if (z > 0.0) U += std::pow(z, degree);
    }
    out[i] = (lower ? L : U) / (L + U);
  }
  return out;
}

NumericVector LPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  if (degree>0) {
    NumericVector lpm_output = LPM_CPv(degree, target, variable);
    NumericVector upm_output = UPM_CPv(degree, target, variable);
    return rescale_unsafe_ratio(lpm_output, upm_output, degree, target, variable, true);
  } else {
    return LPM_CPv(degree, target, variable);
  }
}
NumericVector UPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  if (degree>0) {
    NumericVector lpm_output = LPM_CPv(degree, target, variable);
    NumericVector upm_output = UPM_CPv(degree, target, variable);
    return rescale_unsafe_ratio(lpm_output, upm_output, degree, target, variable, false);
  } else {
    return UPM_CPv(degree, target, variable);
  }
}

double CoUPM_C(
    const double &degree_x, const double &degree_y,
    const RVector<double> &x, const RVector<double> &y,
    const double &target_x, const double &target_y
){
  size_t n_x = x.size(), n_y = y.size();
  size_t max_size = (n_x>n_y ? n_x : n_y);
  size_t min_size = (n_x<n_y ? n_x : n_y);
  if (min_size<=0)
    return 0;
  
  double out=0;
  bool d_x_0=(degree_x==0);
  bool d_y_0=(degree_y==0);
  bool x_is_int=nns_is_integer_degree(degree_x);
  bool y_is_int=nns_is_integer_degree(degree_y);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    double y1=(y[i]-target_y);
    
    if(d_x_0) x1 = (x1 > 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_y_0) y1 = (y1 > 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(!d_x_0){
      if(x_is_int) x1 = nns_int_pow(x1, static_cast<int>(degree_x));
      else x1 = std::pow(x1, degree_x);
    }
    if(!d_y_0){
      if(y_is_int) y1 = nns_int_pow(y1, static_cast<int>(degree_y));
      else y1 = std::pow(y1, degree_y);
    }
    out += x1 * y1;
  }
  return out/max_size;
}

double CoLPM_C(
    const double &degree_x, const double &degree_y,
    const RVector<double> &x, const RVector<double> &y,
    const double &target_x, const double &target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (min_size<=0)
    return 0;
  double out=0;
  bool d_x_0=(degree_x==0);
  bool d_y_0=(degree_y==0);
  bool x_is_int=nns_is_integer_degree(degree_x);
  bool y_is_int=nns_is_integer_degree(degree_y);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    double y1=(target_y-y[i]);
    
    if(d_x_0) x1 = (x1 >= 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_y_0) y1 = (y1 >= 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(!d_x_0){
      if(x_is_int) x1 = nns_int_pow(x1, static_cast<int>(degree_x));
      else x1 = std::pow(x1, degree_x);
    }
    if(!d_y_0){
      if(y_is_int) y1 = nns_int_pow(y1, static_cast<int>(degree_y));
      else y1 = std::pow(y1, degree_y);
    }
    out += x1 * y1;
  }
  return out/max_size;
}

double DLPM_C(
    const double &degree_lpm, const double &degree_upm,
    const RVector<double> &x, const RVector<double> &y,
    const double &target_x, const double &target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (min_size<=0)
    return 0;
  double out=0;
  bool dont_use_pow_lpm=nns_is_integer_degree(degree_lpm),
    dont_use_pow_upm=nns_is_integer_degree(degree_upm),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    double y1=(target_y-y[i]);
    
    if(d_upm_0) x1 = (x1 > 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_lpm_0) y1 = (y1 >= 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(dont_use_pow_lpm && dont_use_pow_upm){
      if(!d_upm_0) x1 = nns_int_pow(x1, static_cast<int>(degree_upm));
      if(!d_lpm_0) y1 = nns_int_pow(y1, static_cast<int>(degree_lpm));
      out += x1 * y1;
    } else if(dont_use_pow_lpm && !dont_use_pow_upm){
      if(!d_lpm_0) y1 = nns_int_pow(y1, static_cast<int>(degree_lpm));
      out += std::pow(x1, degree_upm) * y1;
    } else if(dont_use_pow_upm && !dont_use_pow_lpm){
      if(!d_upm_0) x1 = nns_int_pow(x1, static_cast<int>(degree_upm));
      out += x1 * std::pow(y1, degree_lpm);
    } else out += std::pow(x1, degree_upm) * std::pow(y1, degree_lpm);
  }
  return out/max_size;
}

double DUPM_C(
    const double &degree_lpm, const double &degree_upm,
    const RVector<double> &x, const RVector<double> &y,
    const double &target_x, const double &target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (min_size<=0)
    return 0;
  double out=0;
  
  bool dont_use_pow_lpm=(nns_is_integer_degree(degree_lpm)),
    dont_use_pow_upm=(nns_is_integer_degree(degree_upm)),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    double y1=(y[i]-target_y);
    
    if(d_lpm_0) x1 = (x1 >= 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_upm_0) y1 = (y1 > 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(dont_use_pow_lpm && dont_use_pow_upm){
      if(!d_lpm_0) x1 = nns_int_pow(x1, static_cast<int>(degree_lpm));
      if(!d_upm_0) y1 = nns_int_pow(y1, static_cast<int>(degree_upm));
      out += x1 * y1;
    } else if(dont_use_pow_lpm && !dont_use_pow_upm){
      if(!d_lpm_0) x1 = nns_int_pow(x1, static_cast<int>(degree_lpm));
      out += x1 * std::pow(y1, degree_upm);
    } else if(dont_use_pow_upm && !dont_use_pow_lpm){
      if(!d_upm_0) y1 = nns_int_pow(y1, static_cast<int>(degree_upm));
      out += std::pow(x1, degree_lpm) * y1;
    } else out += std::pow(x1, degree_lpm) * std::pow(y1, degree_upm);
  }
  return out/max_size;
}

#define NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS, LPM_DEGREE_VARIABLE, UPM_DEGREE_VARIABLE) \
size_t target_x_size=target_x.size();                                                               \
size_t target_y_size=target_y.size();                                                               \
size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);                   \
NumericVector output = NumericVector(max_target_size);                                              \
WORKER_CLASS tmp_func(LPM_DEGREE_VARIABLE, UPM_DEGREE_VARIABLE, x, y, target_x, target_y, output);  \
nns_parallel_for(0, output.size(), tmp_func, double(output.size()) * double(x.size()));             \
return(output);

NumericVector CoLPM_CPv(
    const double &degree_x, const double &degree_y,
    const NumericVector &x, const NumericVector &y,
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoLPM_Worker, degree_x, degree_y);
}
NumericVector CoUPM_CPv(
    const double &degree_x, const double &degree_y,
    const NumericVector &x, const NumericVector &y,
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoUPM_Worker, degree_x, degree_y);
}
NumericVector DLPM_CPv(
    const double &degree_lpm, const double &degree_upm,
    const NumericVector &x, const NumericVector &y,
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(DLPM_Worker, degree_lpm, degree_upm);
}
NumericVector DUPM_CPv(
    const double &degree_lpm, const double &degree_upm,
    const NumericVector &x, const NumericVector &y,
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(DUPM_Worker, degree_lpm, degree_upm);
}

// Retained for absolute backward compatibility with internal single-pair calls
void PMMatrix_Cv(
    const double &degree_lpm,
    const double &degree_upm,
    const RMatrix<double>::Column &x,
    const RMatrix<double>::Column &y,
    const double &target_x,
    const double &target_y,
    const bool &pop_adj,
    const double &adjust,
    const size_t &rows,
    double &coLpm,
    double &coUpm,
    double &dLpm,
    double &dUpm,
    double &covMat
){
  RVector<double> x_rvec(x);
  RVector<double> y_rvec(y);
  
  coLpm = 0.0;
  coUpm = 0.0;
  dLpm = 0.0;
  dUpm = 0.0;
  covMat=0;
  if(rows == 0)
    return;
  
  bool lpm_is_int = nns_is_integer_degree(degree_lpm);
  bool upm_is_int = nns_is_integer_degree(degree_upm);
  for(size_t i=0; i<rows; i++){
    double x_lower = lower_component(target_x - x_rvec[i], degree_lpm, lpm_is_int);
    double x_upper = upper_component(x_rvec[i] - target_x, degree_upm, upm_is_int);
    double y_lower = lower_component(target_y - y_rvec[i], degree_lpm, lpm_is_int);
    double y_upper = upper_component(y_rvec[i] - target_y, degree_upm, upm_is_int);
    
    coLpm += x_lower * y_lower;
    coUpm += x_upper * y_upper;
    dLpm += x_upper * y_lower;
    dUpm += x_lower * y_upper;
  }
  
  double inv_rows = 1.0 / static_cast<double>(rows);
  coLpm *= inv_rows;
  coUpm *= inv_rows;
  dLpm *= inv_rows;
  dUpm *= inv_rows;
  
  if(pop_adj && rows > 1 && degree_lpm > 0 && degree_upm > 0){
    coLpm *= adjust;
    coUpm *= adjust;
    dLpm *= adjust;
    dUpm *= adjust;
  }
  covMat = coUpm + coLpm - dUpm - dLpm;
}

// ============================================================================
// ULTRA-OPTIMIZED TENSORIZED MULTIVARIATE INTERNALS
// ============================================================================

// Worker 1: Compute Deviation Matrices exactly ONCE per element.
// Perfectly column-contiguous, cache-friendly SIMD streaming.
struct PrecomputeDeviationsWorker : public Worker {
  const RMatrix<double> variable;
  const RVector<double> target;
  double degree_lpm;
  double degree_upm;
  bool lpm_is_int;
  bool upm_is_int;
  
  RMatrix<double> D_lower;
  RMatrix<double> D_upper;
  
  PrecomputeDeviationsWorker(const NumericMatrix& variable_, const NumericVector& target_,
                             double degree_lpm_, double degree_upm_,
                             NumericMatrix& D_lower_, NumericMatrix& D_upper_)
    : variable(variable_), target(target_),
      degree_lpm(degree_lpm_), degree_upm(degree_upm_),
      lpm_is_int(nns_is_integer_degree(degree_lpm_)), upm_is_int(nns_is_integer_degree(degree_upm_)),
      D_lower(D_lower_), D_upper(D_upper_) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    size_t rows = variable.nrow();
    for (std::size_t j = begin; j < end; ++j) {
      double t_j = target[j];
      for (size_t i = 0; i < rows; ++i) {
        double val = variable(i, j);
        D_lower(i, j) = lower_component(t_j - val, degree_lpm, lpm_is_int);
        D_upper(i, j) = upper_component(val - t_j, degree_upm, upm_is_int);
      }
    }
  }
};

// Worker 2: Blistering Fused Matrix Multiplication (t(D) %*% D)
// Completely stripped of all conditions, branching, and pow() calls.
struct FusedMatrixMultiplicationWorker : public Worker {
  const RMatrix<double> D_lower;
  const RMatrix<double> D_upper;
  bool apply_adj;
  double adjust;
  size_t rows;
  
  RMatrix<double> coLpm;
  RMatrix<double> coUpm;
  RMatrix<double> dLpm;
  RMatrix<double> dUpm;
  RMatrix<double> covMat;
  
  FusedMatrixMultiplicationWorker(const NumericMatrix& D_lower_, const NumericMatrix& D_upper_,
                                  bool apply_adj_, double adjust_, size_t rows_,
                                  NumericMatrix& coLpm_, NumericMatrix& coUpm_,
                                  NumericMatrix& dLpm_, NumericMatrix& dUpm_, NumericMatrix& covMat_)
    : D_lower(D_lower_), D_upper(D_upper_), apply_adj(apply_adj_), adjust(adjust_), rows(rows_),
      coLpm(coLpm_), coUpm(coUpm_), dLpm(dLpm_), dUpm(dUpm_), covMat(covMat_) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    size_t cols = D_lower.ncol();
    double inv_rows = 1.0 / static_cast<double>(rows);
    
    for (std::size_t i = begin; i < end; ++i) {
      // PM.matrix quadrant symmetry:
      //   CUPM(i,j) = CUPM(j,i)
      //   CLPM(i,j) = CLPM(j,i)
      //   DUPM(i,j) = DLPM(j,i)
      //   DLPM(i,j) = DUPM(j,i)
      // Therefore compute only the upper triangle and cross-mirror DUPM/DLPM.
      for (std::size_t j = i; j < cols; ++j) {
        double sum_cupm = 0.0;
        double sum_clpm = 0.0;
        double sum_dupm = 0.0;
        double sum_dlpm = 0.0;
        
        // Loop fusion: Compute all 4 co-moment quadrants in a single hot-cache row scan.
        for (size_t k = 0; k < rows; ++k) {
          double u_i = D_upper(k, i);
          double l_i = D_lower(k, i);
          double u_j = D_upper(k, j);
          double l_j = D_lower(k, j);
          
          sum_cupm += u_i * u_j;
          sum_clpm += l_i * l_j;
          sum_dupm += l_i * u_j;
          sum_dlpm += u_i * l_j;
        }
        
        sum_cupm *= inv_rows;
        sum_clpm *= inv_rows;
        sum_dupm *= inv_rows;
        sum_dlpm *= inv_rows;
        
        if (apply_adj) {
          sum_cupm *= adjust;
          sum_clpm *= adjust;
          sum_dupm *= adjust;
          sum_dlpm *= adjust;
        }
        
        double cov_ij = sum_cupm + sum_clpm - sum_dupm - sum_dlpm;
        
        coUpm(i, j) = sum_cupm;
        coLpm(i, j) = sum_clpm;
        dUpm(i, j)  = sum_dupm;
        dLpm(i, j)  = sum_dlpm;
        covMat(i, j) = cov_ij;
        
        if (j != i) {
          coUpm(j, i) = sum_cupm;
          coLpm(j, i) = sum_clpm;
          
          // Crossed mirror, not ordinary symmetry.
          dUpm(j, i)  = sum_dlpm;
          dLpm(j, i)  = sum_dupm;
          covMat(j, i) = cov_ij;
        }
      }
    }
  }
};

// Normalised PM.matrix cell (i, j) for a pair whose plain quadrant total is
// unsafe (zero, under- or overflowed). out = {cupm, dupm, dlpm, clpm}.
//
// With LPM_degree > 0 a tie (x == t) contributes 0 to every quadrant, so a row
// contributes nothing iff it ties the target in column i or j. If every row
// does, the total is exactly 0 (0/0): return the degree-0 cells (lower side
// inclusive, upper strict), as LPM.ratio does for a constant variable and the
// nD functions do for a zero total; they partition the sample, so the cell
// still sums to 1. With LPM_degree = 0 a tie counts as lower and the total is
// never exactly 0.
//
// Otherwise the total under- or overflowed. Column c's deviations are divided by
// s_c = max_k |x_kc - t_c|, and quadrant q, which raises column i to a_q and
// column j to b_q, is multiplied back by 2^(a_q log2 s_i + b_q log2 s_j)
// relative to the largest such factor. With equal degrees the four factors are
// the same expression, so they cancel exactly and the cell is the plain
// per-column rescaled ratio. If the rescaled sums still under- or overflow
// (very high degree), the terms are summed in the log domain.
//
// Accuracy: relative error of order degree * eps, not exact. The exponent
// a_q log2 s_i + b_q log2 s_j (about 1e3 to 1e4 in magnitude at degree 1000)
// is rounded before exp2, and the rounding of |d| / s is raised to the degree.
// Splitting s into mantissa and integer exponent would remove the first source
// but not the second. Both vanish only for dyadic data (s a power of 2 and
// |d| / s exact).
static void PMMatrix_normalized_rescue(const NumericMatrix &variable,
                                       const NumericVector &target,
                                       const std::size_t i, const std::size_t j,
                                       const double degree_lpm, const double degree_upm,
                                       const std::vector<double> &scale,
                                       double out[4]) {
  const std::size_t n = variable.nrow();
  const bool lpm_is_int = nns_is_integer_degree(degree_lpm);
  const bool upm_is_int = nns_is_integer_degree(degree_upm);
  const double ti = target[i], tj = target[j];
  const double si = scale[i], sj = scale[j];
  // log2 of a zero scale is never used with a nonzero exponent: an all-tied
  // column has z = 0, which only contributes through a degree-0 indicator.
  const double lsi = si > 0.0 ? std::log2(si) : 0.0;
  const double lsj = sj > 0.0 ? std::log2(sj) : 0.0;
  // Quadrant exponents (column i, column j): cupm (U, U), dupm (L, U),
  // dlpm (U, L), clpm (L, L).
  const double a[4] = {degree_upm, degree_lpm, degree_upm, degree_lpm};
  const double b[4] = {degree_upm, degree_upm, degree_lpm, degree_lpm};
  double offset[4];
  for (int q = 0; q < 4; ++q)
    offset[q] = (a[q] == 0.0 ? 0.0 : a[q] * lsi) + (b[q] == 0.0 ? 0.0 : b[q] * lsj);

  double c0[4] = {0.0, 0.0, 0.0, 0.0};
  double S[4] = {0.0, 0.0, 0.0, 0.0};
  bool any_term = false;
  for (std::size_t k = 0; k < n; ++k) {
    const double di = variable(k, i) - ti, dj = variable(k, j) - tj;
    const bool up_i = di > 0.0, up_j = dj > 0.0;
    c0[up_i ? (up_j ? 0 : 2) : (up_j ? 1 : 3)] += 1.0;
    const double zi = di == 0.0 ? 0.0 : std::fabs(di) / si;
    const double zj = dj == 0.0 ? 0.0 : std::fabs(dj) / sj;
    const double li = lower_component(up_i ? -zi : zi, degree_lpm, lpm_is_int);
    const double ui = upper_component(up_i ? zi : -zi, degree_upm, upm_is_int);
    const double lj = lower_component(up_j ? -zj : zj, degree_lpm, lpm_is_int);
    const double uj = upper_component(up_j ? zj : -zj, degree_upm, upm_is_int);
    const double term[4] = {ui * uj, li * uj, ui * lj, li * lj};
    for (int q = 0; q < 4; ++q) S[q] += term[q];
    if ((up_i || degree_lpm == 0.0 || di != 0.0) && (up_j || degree_lpm == 0.0 || dj != 0.0))
      any_term = true;
  }

  if (!any_term) {
    for (int q = 0; q < 4; ++q) out[q] = c0[q] / double(n);
    return;
  }

  double omax = -std::numeric_limits<double>::infinity();
  for (int q = 0; q < 4; ++q) if (S[q] > 0.0) omax = std::max(omax, offset[q]);
  double w[4], total = 0.0;
  for (int q = 0; q < 4; ++q) {
    w[q] = S[q] > 0.0 ? S[q] * std::exp2(offset[q] - omax) : 0.0;
    total += w[q];
  }

  if (nns_unsafe_total(total)) {
    // Log domain: log2 of each nonzero term, shifted by the largest. log2|d|
    // - log2 s rather than log2(|d| / s): the quotient underflows to 0 when a
    // column's deviations span more than 2^1074 (a denormal next to 1e300).
    std::vector<double> lg(4 * n, -std::numeric_limits<double>::infinity());
    double lmax = -std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < n; ++k) {
      const double di = variable(k, i) - ti, dj = variable(k, j) - tj;
      const bool up_i = di > 0.0, up_j = dj > 0.0;
      // Side of each column for quadrant q: true = upper.
      const bool side_i[4] = {true, false, true, false};
      const bool side_j[4] = {true, true, false, false};
      for (int q = 0; q < 4; ++q) {
        if (side_i[q] != up_i || side_j[q] != up_j) continue;
        double l = offset[q];
        bool zero = false;
        if (a[q] != 0.0) { if (di == 0.0) zero = true; else l += a[q] * (std::log2(std::fabs(di)) - lsi); }
        if (b[q] != 0.0) { if (dj == 0.0) zero = true; else l += b[q] * (std::log2(std::fabs(dj)) - lsj); }
        if (zero) continue;
        lg[4 * k + q] = l;
        lmax = std::max(lmax, l);
      }
    }
    total = 0.0;
    for (int q = 0; q < 4; ++q) w[q] = 0.0;
    for (std::size_t k = 0; k < n; ++k)
      for (int q = 0; q < 4; ++q)
        if (lg[4 * k + q] > -std::numeric_limits<double>::infinity())
          w[q] += std::exp2(lg[4 * k + q] - lmax);
    for (int q = 0; q < 4; ++q) total += w[q];
  }
  for (int q = 0; q < 4; ++q) out[q] = w[q] / total;
}

// [[Rcpp::export]]
List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const NumericVector &target,
    const NumericMatrix &variable,
    const bool &pop_adj,
    const bool &norm
) {
  nns_check_degree(LPM_degree, "LPM_degree");
  nns_check_degree(UPM_degree, "UPM_degree");
  size_t variable_cols = variable.cols();
  size_t target_length = target.size();
  if(variable_cols != target_length){
    Rcpp::stop("variable matrix cols != target vector length");
    return List::create();
  }
  
  size_t rows = variable.rows();
  if (rows == 0) Rcpp::stop("`variable` must have at least one row");
  
  // 1. Allocate continuous intermediate deviation matrices
  NumericMatrix D_lower(rows, variable_cols);
  NumericMatrix D_upper(rows, variable_cols);
  
  // 2. Step 1: Precompute all element deviation components in parallel
  PrecomputeDeviationsWorker precalc_engine(variable, target, LPM_degree, UPM_degree, D_lower, D_upper);
  nns_parallel_for(0, variable_cols, precalc_engine, double(rows) * double(variable_cols));
  
  // 3. Allocate final return matrix structures
  NumericMatrix coLpm(variable_cols, variable_cols);
  NumericMatrix coUpm(variable_cols, variable_cols);
  NumericMatrix dLpm(variable_cols, variable_cols);
  NumericMatrix dUpm(variable_cols, variable_cols);
  NumericMatrix covMat(variable_cols, variable_cols);
  
  // 4. Determine population adjustment configurations
  double adjust = 1.0;
  if (pop_adj && rows > 1) {
    adjust = static_cast<double>(rows) / static_cast<double>(rows - 1);
  }
  bool apply_adj = pop_adj && rows > 1 && LPM_degree > 0 && UPM_degree > 0;
  
  // 5. Step 2: High-speed matrix contraction loops across available cores
  FusedMatrixMultiplicationWorker matrix_engine(D_lower, D_upper, apply_adj, adjust, rows,
                                                coLpm, coUpm, dLpm, dUpm, covMat);
  nns_parallel_for(0, variable_cols, matrix_engine, double(rows) * double(variable_cols) * double(variable_cols));
  
  // 6. Apply cellular normalization adjustments if requested.
  // Preserve the same cross-transpose relationship for DUPM and DLPM.
  // A cell with an unsafe total (zero, under- or overflowed) is recomputed by
  // PMMatrix_normalized_rescue; s_c = max_k |x_kc - t_c| is computed on demand.
  if (norm) {
    std::vector<double> scale;
    for (size_t i = 0; i < variable_cols; ++i) {
      for (size_t j = i; j < variable_cols; ++j) {
        double cupm_ij = coUpm(i, j);
        double dupm_ij = dUpm(i, j);
        double dlpm_ij = dLpm(i, j);
        double clpm_ij = coLpm(i, j);
        double total = cupm_ij + dupm_ij + dlpm_ij + clpm_ij;
        
        if (!nns_unsafe_total(total)) {
          cupm_ij /= total;
          dupm_ij /= total;
          dlpm_ij /= total;
          clpm_ij /= total;
        } else {
          if (scale.empty()) {
            scale.assign(variable_cols, 0.0);
            for (size_t c = 0; c < variable_cols; ++c)
              for (size_t k = 0; k < rows; ++k)
                scale[c] = std::max(scale[c], std::fabs(variable(k, c) - target[c]));
          }
          double cell[4];
          PMMatrix_normalized_rescue(variable, target, i, j, LPM_degree, UPM_degree, scale, cell);
          cupm_ij = cell[0];
          dupm_ij = cell[1];
          dlpm_ij = cell[2];
          clpm_ij = cell[3];
        }
        
        double cov_ij = cupm_ij + clpm_ij - dupm_ij - dlpm_ij;
        
        coUpm(i, j) = cupm_ij;
        coLpm(i, j) = clpm_ij;
        dUpm(i, j)  = dupm_ij;
        dLpm(i, j)  = dlpm_ij;
        covMat(i, j) = cov_ij;
        
        if (j != i) {
          coUpm(j, i) = cupm_ij;
          coLpm(j, i) = clpm_ij;
          
          // Crossed mirror after normalization too.
          dUpm(j, i)  = dlpm_ij;
          dLpm(j, i)  = dupm_ij;
          covMat(j, i) = cov_ij;
        }
      }
    }
  }
  
  // 7. Shape attribute text allocations
  rownames(coLpm) = colnames(variable);
  colnames(coLpm) = colnames(variable);
  
  rownames(coUpm) = colnames(variable);
  colnames(coUpm) = colnames(variable);
  
  rownames(dLpm) = colnames(variable);
  colnames(dLpm) = colnames(variable);
  
  rownames(dUpm) = colnames(variable);
  colnames(dUpm) = colnames(variable);
  
  rownames(covMat) = colnames(variable);
  colnames(covMat) = colnames(variable);
  
  return(
    List::create(
      Named("cupm") = coUpm,
      Named("dupm") = dUpm,
      Named("dlpm") = dLpm,
      Named("clpm") = coLpm,
      Named("cov.matrix") = covMat
    )
  );
}

# NNS 14.0 partial-moment audit (unreleased)

* Partial-moment functions (`LPM`, `UPM`, their ratios, `Co.*`, `D.*`,
  `PM.matrix`, the `_nD` functions and `NNS.CDF`) now reject an infinite
  degree with "`degree` must be finite." (named after the argument) instead
  of returning Inf, 0 or NaN.
* `LPM()`/`UPM()`/`LPM.ratio()`/`UPM.ratio()` with more than 32 targets, and
  univariate `NNS.CDF()` with more than 32 observations, no longer round a
  degree within 1e-12 of an integer to that integer, so they agree with the
  32-or-fewer-target path (e.g. `UPM(2 + 1e-13, ...)` is no longer exactly the
  degree-2 value). For `LPM()` and the ratios, a tiny positive degree such as
  `1e-13` now excludes ties at the target, as it already did with 32 or fewer
  targets, instead of returning the degree-0 value.
* Integer degrees above 32 are evaluated with `pow()` rather than a product of
  `degree` factors, so a huge integer degree no longer takes time proportional
  to the degree. Results for degrees up to 32 are unchanged; above 32 they can
  differ in the last bits.
* `LPM()`, `UPM()`, `LPM.ratio()`, `UPM.ratio()` with more than 32 targets,
  and univariate `NNS.CDF()` with more than 32 observations, return `Inf` (or
  the corresponding 0/1 ratio) instead of `NaN` when an intermediate power sum
  overflows, matching the 32-or-fewer-target path. (The ratios no longer
  return `NaN` from `Inf/Inf` either; see the scale-free ratio entry below.)
* `LPM()` and `UPM()` require `excess_ret` to be a single `TRUE` or `FALSE`
  (numeric 0/1 still work); `NA`, vectors and strings now give "`excess_ret`
  must be TRUE or FALSE." for every target length. `NA` used to run the
  excess-return path for a scalar target.
* `LPM()` and `UPM()` with `target = NULL` now use `mean(variable)`, like
  `LPM.ratio()`, `Co.LPM()` and `PM.matrix()`, instead of returning
  `numeric(0)`. A zero-length target still returns `numeric(0)` (with
  `excess_ret = FALSE`).
* Factors are rejected with "`x` must be numeric, not a factor." (named after
  the argument) instead of being computed on their level codes, as data or
  target, in `LPM()`, `UPM()`, `LPM.ratio()`, `UPM.ratio()`, `Co.LPM()`,
  `Co.UPM()`, `D.LPM()`, `D.UPM()`, the `_nD` functions, `NNS.CDF()`,
  `LPM.VaR()`, `UPM.VaR()`, `NNS.moments()` and the stochastic dominance
  functions (`NNS.FSD()`, `NNS.SSD()`, `NNS.TSD()` and their `.uni`
  versions), also when the factor is the column of a one-column data frame.
  A factor column of multi-column data is rejected the same way in the `_nD`
  functions, multivariate `NNS.CDF()` and `LPM.VaR()`/`UPM.VaR()`, where it
  used to fail with "`data` must be finite." or be flattened to level codes.
* A target that is not `NULL`, numeric or logical (a string, list, data frame
  or factor) is now an error in `LPM.ratio()`, `UPM.ratio()`, `Co.LPM()`,
  `Co.UPM()`, `D.LPM()`, `D.UPM()` and `PM.matrix()` instead of silently
  becoming `mean(variable)`. A logical target is coerced as `LPM()` does
  (`TRUE` is 1, `NA` gives "`target` must be finite."). `PM.matrix()` keeps
  accepting exactly `target = "mean"` and `"means"` as aliases for `NULL`.
* `LPM()`, `UPM()`, `LPM.VaR()`, `UPM.VaR()`, `NNS.moments()`, `NNS.CDF()`
  targets and the stochastic dominance functions now reject character and
  list data ("`x` must be numeric.") instead of parsing or coercing it.
  For `LPM()`/`UPM()` this includes a `NULL` `variable`, which used to give
  `NaN` and is now an error as in `LPM.ratio()`; `numeric(0)` still gives
  `NaN`.
* `LPM.VaR()` and `UPM.VaR()` flatten any data frame, as they already did a
  tibble or data.table, so a plain `data.frame` is accepted instead of
  being an error.
* `LPM.VaR()` and `UPM.VaR()` reject infinite values in `x` with "`x` must be
  finite." at every degree. Degrees 1-4 used to fail with "missing value
  where TRUE/FALSE needed", and degree 0 returned a quantile. `NA`/`NaN` are
  still dropped.
* `LPM.VaR()` and `UPM.VaR()` reject a missing `percentile` with
  "`percentile` must not be NA." and a character one with "`percentile` must
  be numeric." at every degree. A missing percentile gave `NA` at degree 0,
  silently 0 at degrees 1-4 when other percentiles were present, and "missing
  value where TRUE/FALSE needed" otherwise. At degrees 1-4, data whose spread
  from the median overflows (such as `c(-1e308, 1e308)`) is centred at the
  midrange instead of failing with the same message.
* `NNS.CDF()` accepts a one-column data frame as univariate data and a vector
  of targets for univariate data (evaluated element-wise); a missing target
  gives "`target` must be finite.".
* `PM.matrix()` with a zero-row `variable` errors with "`variable` must have
  at least one row" for every target, instead of returning `list()` or a
  non-finite target error; `Co.LPM_nD()`, `Co.UPM_nD()` and `DPM_nD()` on
  zero-row data error with "`data` must have at least one row", as the batch
  backend used by multivariate `NNS.CDF()` already did, instead of returning
  `NaN` or 0.
* `LPM()`/`UPM()` with `excess_ret = TRUE` index `variable` with `R_xlen_t`,
  so vectors longer than `.Machine$integer.max` are no longer truncated.
* `LPM.ratio()` and `UPM.ratio()` are scale-free at any degree: where the raw
  sums underflow or overflow (high degree, or deviations of extreme
  magnitude) the ratio is recomputed on deviations divided by
  `max(abs(variable - target))` instead of returning `NaN` (e.g.
  `LPM.ratio(400, 0, c(-0.01, 0.01))` is 0.5). This also removes the `NaN`
  curves of univariate `NNS.CDF()` at high degree and the wrong roots and
  "NA/NaN replaced" warnings of `LPM.VaR()`/`UPM.VaR()` for non-integer or
  > 4 degrees. Results with a safe total are unchanged.
* `LPM.VaR()`/`UPM.VaR()` with integer degree 1 to 4 solve on data rescaled by
  an exact power of two when the deviations are below 2^-200 or above 2^200,
  so data of scale 1e-200 no longer returns the same value for every
  percentile. A degree above `.Machine$integer.max` no longer warns "NAs
  introduced by coercion to integer range".
* `NNS.moments()` computes skewness and kurtosis from deviations rescaled by an
  exact power of two at extreme scales, so data of scale 1e-200 or 1e200 gives
  finite values instead of `NaN` (the variance is still `mean((x - mean(x))^2)`,
  which can underflow to 0).
* Normalised `Co.LPM_nD()`, `Co.UPM_nD()` and `DPM_nD()` at degree > 0, and
  multivariate `NNS.CDF()` at degree > 0, no longer return 0 for all three
  cells when the total is 0. If every observation ties the target in some
  coordinate (e.g. a constant column) the result is the degree-0 value, as
  `LPM.ratio()` does for a constant variable, so the cells sum to 1 and the
  degree-1 CDF of data with a constant column is no longer 0 everywhere. A
  total that underflowed or overflowed is recomputed on per-column rescaled
  deviations.
* Multivariate `NNS.CDF(type = "survival")` at degree > 0 is computed by the
  batched `Co.UPM_nD()` directly (bit-identical on ordinary data), so its
  zero-total fallback is the strict all-above count of `Co.UPM_nD()`.
* `LPM.VaR()`/`UPM.VaR()` at integer degrees 1 to 4 return the root of the
  ratio equation, not a nearby order statistic. The solver used to accept a
  break whenever its ratio was within about 1.8e-12 of the percentile, which
  is not a tolerance on the threshold where the ratio is flat:
  `LPM.VaR(1e-12, 4, c(0, 1))` was 0 and `UPM.VaR(1e-12, 4, c(0, 1))` was 1
  (root about 0.001 from each end), and `LPM.VaR(0.5 + 1e-13, 2, c(0, 1, 2))`
  was 1. Its prefix power sums were expanded around the median and cancelled
  from degree 3 up for a cluster far from it (a percentile of 1e-21 in a
  cluster at -1e5 returned the wrong order statistic, 1.2 away). The gap is
  now located by the sign of `(1 - p) LPM - p UPM` at the order statistics,
  using sums anchored at the adjacent order statistics (no cancellation,
  accumulated in long double), and solved there from the nearer end with
  `uniroot()`. A root more than 30 binades closer to that end than half the
  gap (a tail percentile) is bracketed by a bound on its offset, so a
  percentile of 1e-300 takes tens of `uniroot()` iterations, not about a
  thousand. `UPM.VaR()` keeps its percentile exact instead of solving at
  `1 - percentile`. On ordinary data the results change only by rounding
  (less than 2e-14 times the largest threshold in the audit's samples), and
  the help no longer states a ratio tolerance. Degree 0 and the
  non-integer-degree solver are unchanged. `LPM.VaR(LPM.ratio(d, x_k, x), d,
  x)` is no longer bit-identical to `x_k` at the extreme order statistics,
  where rounding of the ratio near 0 or 1 moves the root (by about 1e-8 at
  `d = 4` for the second-largest of 300 standard normal draws); it is the
  root of the rounded ratio.
* An invalid `type` in `NNS.FSD()`, `NNS.FSD.uni()`, `NNS.SD.efficient.set()`
  and `NNS.SD.cluster()` (anything but `"discrete"` or `"continuous"`, in any
  letter case) is now an error, "`type` must be \"discrete\" or
  \"continuous\".", at every degree. It used to warn and then fall back to
  the Omega order in `NNS.FSD()` but to discrete FSD everywhere else.
* The stochastic dominance functions reject infinite observations with
  "`x` must be finite." (naming `y` or `data` when that is the argument).
  `NNS.SSD.uni()`, `NNS.TSD.uni()`, continuous `NNS.FSD.uni()` and the SSD/TSD
  efficient set and clusters used to report no dominance whenever a value was
  infinite (a total-loss log return made the asset SSD-efficient although it
  was FSD-dominated); discrete `NNS.FSD.uni()` accepted it. `NNS.FSD()`,
  `NNS.SSD()` and `NNS.TSD()` already failed, but with "`variable` must be
  finite.". `NA`/`NaN` keep the "missing values" error.
* `NNS.FSD()`, `NNS.SSD()` and `NNS.TSD()` compute their plot curves only when
  `plot = TRUE`; the verdict never used them.
* The `NNS.TSD()` plot draws the LPM2 curves against the threshold values, like
  `NNS.SSD()` and `NNS.FSD()`, instead of against their index.
* `PM.matrix(norm = TRUE)` with `LPM_degree > 0` no longer sets all four
  quadrants of a cell to 0 when its total is 0, i.e. when every row ties the
  target in column i or j (a column constant at its target, or `target =
  NULL` with a constant column). The cell takes the `LPM_degree = UPM_degree =
  0` values, ties counting as lower, as `LPM.ratio()` does for a constant
  variable and `Co.LPM_nD()`/`Co.UPM_nD()`/`DPM_nD()` do for a zero total, so
  the four quadrants still sum to 1: such a column has `$cov.matrix` diagonal 1
  (was 0), and its off-diagonal entries are `2 * F_j(t_j) - 1` (was 0), with
  `F_j(t_j)` the share of the other column at or below its target. A cell
  total that underflowed or overflowed (deviations of scale 1e-200, or high
  degree) is recomputed on per-column rescaled deviations, in the log domain
  if needed, instead of returning 0 for every quadrant; with equal degrees the
  result does not depend on the units of the data. Cells with a safe total
  are unchanged bit for bit. `NNS.copula()` drops constant columns before
  calling `PM.matrix()`, so it changes only for a pair of columns in which
  every row ties a user-supplied target in one of the two.

# NNS 13.3 (development changes)

## Statistical targets and interpretation

* Preserve the distinction between frequency probabilities (degree zero),
  directional deviation-area shares (positive degrees), and their induced laws.
  Correct one-sided probability bounds, divergent-moment labels, finite-moment
  conditions, prediction-interval claims, and the scope of dependence/causation.
* ANOVA retains its separate mean, median, and distribution-profile targets.
  Median comparisons now handle CDF jumps, and degenerate area-balance cases
  use explicit agreement semantics. Certainty remains a descriptive score;
  inferential calibration must match the target and sampling design.
* `NNS.FSD(type = "continuous")` keeps the same numerical degree-one criterion
  but labels it Omega/expectile order. Its result strings now use `OMEGA`, not
  `FSD`; callers comparing strings should update. Discrete FSD strings and
  numeric dominance indicators/set/cluster structures are unchanged.

## Regression, classification, and derivatives

* `NNS.reg` and `NNS.stack` accept `type = "REG"` to select numeric
  regression explicitly. Forecasting, derivative, and hazard-smoothing helpers
  use this contract even when a response has only two values; they retain
  fractional predictions on the original scale.
* Nominal classification retains class probabilities in predictor cells and
  combines those probabilities before choosing labels. Numeric class labels are
  not interpolation coordinates. Probability outputs and descriptive prediction
  sets expose the target; factor/character output codes use `class.levels`.
* `NNS.stack` classification blends its two models with weights proportional
  to their squared objective scores (inverse squares when `objective = "min"`),
  as in NNS 13.0, and tunes the binary threshold for that blend. The weight had
  been the value in 0, 0.01, ..., 1 with the best out-of-fold objective, which
  gave one model all the weight when it was ahead by a single out-of-fold row:
  on an iris split with out-of-fold accuracies of 0.942 and 0.949 (130 and 131
  of 138 rows), the weights were 0 and 1, and the stack missed a test row that
  the 13.0 weights (0.496 and 0.504) classify correctly.
* `NNS.reg`, `NNS.stack` and `NNS.boost` match complete predictor
  names in any order and use positional matching when the test and training
  columns share no names. The known Excel CSV UTF-8 first-header artifact
  (`Male` becoming `X...Male`) is repaired before name matching. Arbitrary
  partial renames or misspelled predictors now fail instead of being silently
  treated as positional matches.
* Classification `NNS.stack` reports `NNS.dim.red.threshold` as the
  coefficient cutoff at which `NNS.reg` on the training data keeps the tuned
  number of predictors, as regression stacks do: one number when the class
  indices share it, as two-class indices always do, and otherwise one cutoff
  per class, which `NNS.reg` accepts as `threshold`. It had been `NA`, which
  `NNS.reg` rejected. When all cells are the best neighbour count, the stack
  reports their number rather than `NA`.
* `NNS.reg` treats `threshold = NULL` or `NA` as 0 and `n.best = NA` as the
  default, so `NA` values reported by `NNS.stack` can be passed back.
* Binary numeric coordinates retain distinct terminal cells in both R and
  native regression paths, including dummy columns in mixed predictor designs.
* Stack's multivariate component searches selected predictor coordinates; its
  dimension-reduction component separately uses the synthetic index. Complexity
  selection uses both response orientations without changing the underlying
  directional statistics. These corrections can change fitted predictions.
* Distance kernels use consistent exponential rates and sample rank scales
  across scalar, bulk, and parallel paths. Partial derivatives use admissible
  local support and explicit boundary differences. Unsupported multivariate
  `NNS.CDF(type = "hazard")` requests now fail instead of returning CDF/survival
  as if it were density/survival.
* Stack's multivariate component keeps the same coordinates in every fold, so
  a predictor that is constant in one fold no longer stops the search. A
  declared factor level with no training rows is a valid test value, a single
  test row with a factor predictor is accepted, and categorical `NNS.reg`
  honours an explicit `n.best`.
* Nominal classification reuses each fold's cells and sorted distances for
  every neighbour count and skips fitted classes in point-only calls. It is
  still slower than the former class-code regression: for three classes, two
  predictors and `folds = 2`, `NNS.stack` takes 1.3 s and `NNS.boost` 1.3 s at
  n = 1000 (0.9 s and 0.12 s before), and 5.8 s and 8.0 s at n = 5000
  (2.6 s and 1.2 s before).
* Gravity retains its location blend with corrected quartiles and observed bin
  centers. Mode retains distinct continuous and class behavior. Changes to
  these defaults can affect regression, optimization, and forecasts.

## Forecasting and scenarios

* ARMA validation support is selected from training observations only. Nonlinear
  shrinkage follows the same intended mean blend as linear shrinkage, and a
  single shrunk lag refits recursively. Long-horizon runtime and path bands can
  therefore differ from earlier releases.
* The optimizer's objective describes its returned internal validation forecast.
  `NNS.ARMA.optim` no longer adds the gravity bias shift by default
  (`bias.correction = "none"`). The former rule estimated the shift on the
  validation window and accepted it whenever it improved the objective on that
  same window; it is removed. On rolling-origin forecasts of nine series the
  shift made held-out MASE worse on eight, under the former rule and under
  `bias.correction = "split"`, which estimates it on the first half of the
  window and accepts it on the second and remains an option. Base-model
  tuning still requires independent outer evaluation. Raw/corrected forecasts
  and bias metadata are available for assessment.
* `NNS.ARMA.optim` and the nonlinear `NNS.ARMA` methods run faster. The
  nonlinear step's weighted rise and run were means of a vector about k^3/3
  long for k regression points; they are now computed in closed form, which
  agrees with the former value to rounding. The optimizer computes each
  validation forecast once, partitions no longer build their data frames
  through R, and small partial-moment kernels run on the calling thread; these
  leave results bit-identical. On a synthetic monthly series of 8000
  observations the optimizer took 13 s, against 126 s in 13.3 and 210 s
  before this change; the AirPassengers example takes 0.10 s, against 0.30 s
  and 0.39 s.
* Long series forecast faster with the same results. When every lag is at
  least `h`, no forecast reaches a lag's component series, so `NNS.ARMA`'s
  linear method forecasts several lags directly, as it did one lag, and
  combines them step by step as before; `NNS.ARMA.optim`'s linear search makes
  these calls. `NNS.seas` scores its lags on several threads, and its
  fundamental-period jackknife sorts all components at once. The tides example
  (222,000 six-minute observations, 200 candidate periods) took 11 min on one
  core with 13.0 and 13 min before this change; it now takes 1.0 min, or 0.7
  min with the default cores.
* `NNS.VAR` again completes a variable whose latest values are missing with the
  `NNS.stack` fit on the time index used in NNS 13.0: every missing value takes
  the value of the nearest regression point on the index, the newest
  observation keeps its own regression point, and values after it extend the
  gradient at that boundary. In the 13.1-13.3 betas the multivariate regression
  behind that fit merged the newest observation into the previous cell, stopped
  extending past the data, and later partitioned the index more coarsely, so
  quarterly real GDP in a monthly panel was filled at 21,474 for November 2022
  after an observed 22,273 (13.0: 22,316). On one- and two-month trailing gaps
  in 12 FRED series (origins 2016-2022), the restored fill's median error is
  0.88 mean absolute monthly changes for monthly series and 1.09 for quarterly
  series (13.0: 0.88 and 1.06). When only earlier values are missing they are
  still interpolated linearly.
* VAR keeps its two-stage forecast design and uses prefix-generated predictors
  for chronological validation by default. `validation = "observed"` is the
  compatibility option. Known contemporaneous outcomes and imputed future
  targets are not used as deployment-validation inputs or ground truth.
* `NNS.seas` reports the fundamental period as `best.period`: the smallest
  period that divides the lowest-scoring lag (or nearly divides it, for a
  non-integer period) and scores within two jackknife standard errors of it,
  from paired delete-one differences.
  The lowest score often fell on a multiple of the true period; a monthly sine
  with noise now returns 12 instead of 36 to 108, and `NNS.ARMA`'s automatic
  lag follows. The period is listed first in `all.periods` and `periods`.
  `NNS.ARMA` does not run the permutation test; pass
  `NNS.seas(x, p.value = TRUE)$best.period` as `seasonal.factor` for a
  significance-filtered lag.
* Seasonality adds absolute dispersion for interval-scale data while retaining
  the historical automatic CV/Gini default. Selected predictive periods need not be fundamental
  physical periods; permutation interpretations require an appropriate null.
* Monte Carlo exponent grids stay within requested bounds and reject invalid
  exponents. Constant/singleton superiority groups and exactly zero bootstrap
  residuals are handled explicitly. Bootstrap residual expansion is separately
  configurable. The default forced-mean scale is an iid reference, not a
  universal serial-dependence uncertainty calibration.

## Numerical, documentation, and reproducibility repairs

* Fix mixed-degree co-partial moments, validate paired lengths before native
  parallel work, and make integer-degree VaR root stopping scale invariant.
* Fix target hazard evaluation, stochastic-dominance cluster ordering, residual
  clusters, efficient-set order dependence, and caller column-name mutation.
* Replace unsupported example rankings and inference claims with explicit
  targets and corrected executable workflows. Forecast examples reserve test
  outcomes; fitted-contrast uncertainty includes model refitting. Add revised
  editable companions for the historical PDFs, preserve originals, and record
  that the original editable sources were not recovered.

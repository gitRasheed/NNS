# Exact power-of-two factor for scale-free computations on data of extreme
# magnitude. Returns 1 when `s` (a nonnegative scale such as max|x - centre|) is
# zero, non-finite or inside [lo, hi]; otherwise a power of two k with s * k in
# [0.5, 1). Multiplying by a power of two is exact, and inside the default window
# powers up to degree 4 neither underflow nor overflow, so ordinary data takes
# the k = 1 path and its results are unchanged.
.nns_pow2_rescale <- function(s, lo = 2^-200, hi = 2^200) {
  if (!is.finite(s) || s == 0 || (s >= lo && s <= hi)) return(1)
  2^min(max(-(floor(log2(s)) + 1), -1000), 1000)
}

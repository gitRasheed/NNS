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
  overflows, matching the 32-or-fewer-target path; `Inf/Inf` remains `NaN`.
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

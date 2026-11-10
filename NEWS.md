# BayesChange (development version)

## New features
- New functions `params_uni()`, `params_multi()` and `params_epi()` create the model parameters passed to `detect_cp()` and `clust_cp()`. Their arguments can be auto-completed, have documented defaults and are validated; `print()` shows all values and which ones are defaults. Plain lists are still accepted, and unknown names now raise a warning instead of being silently ignored.
- Model parameters and sampler settings are now separate: `q`, `B`, `L` and `alpha_SM` are only arguments of `detect_cp()`/`clust_cp()`. A `q` given in `params` is still used by `clust_cp()` for univariate time series, with a deprecation warning; other sampler settings in `params` are ignored with a warning.

## Bug fixes
- `user_seed` now fully determines the output of `detect_cp()` and `clust_cp()`. Previously only the GSL generator was seeded, while most draws came from R's random number generator, so two calls with the same `user_seed` could return different results. R's random number generator state is restored on exit; `user_seed = NULL` uses the current state of R's generator.
- Fixed the Metropolis-Hastings update of the autocorrelation parameter `phi` in `detect_cp()` for time series: the multivariate log-gamma function truncated its argument to an integer, the log full conditional evaluated to `-Inf` and every proposal was accepted. The same truncation also affected the integrated likelihood of multivariate time series.
- Fixed the integrated likelihood of univariate time series. It now equals the closed-form marginal of the model in Corradin, Danese & Ongaro (2022) (Eq. (11) with d = 1). Previously the function used by `clust_cp()` and in the split-merge steps had `log(1 - phi)` instead of `log(1 + phi)`, and blocks of length 1 and 2 were handled incorrectly.
- The update of `phi` in `detect_cp()` for univariate time series now uses the prior given by `a`, `b`, `c` (it previously used fixed hyperparameters), and its full conditional no longer omits a `log` in the `k_n` term (this also affected multivariate time series).
- Fixed the importance-sampling estimate of the normalisation constant used by `clust_cp()` (all kernels): the estimate was biased upwards by about `2 * log(B)` and used an incorrect proposal density, which made the algorithm favour too many clusters.
- Fixed a crash (`Mat::operator(): index out of bounds`) in `clust_cp()` for time series, occurring when the random starting partition had as many groups as observations; the starting partition no longer has one group too many.
- The acceleration step of `clust_cp()` for time series now updates the order of the correct cluster (it was indexed by observation instead of by cluster).
- `clust_cp()` for univariate time series now uses its argument `q` (it was ignored in favour of `params$q`, which is still accepted).
- `summary()` for `DetectCpObj` objects reported the average number of blocks (change points + 1) as the average number of change points.
- Computational time is now recorded also when `print_progress = FALSE` (it was reported as 0).
- Random starting orders no longer contain empty blocks.
- `posterior_estimate()` for `DetectCpObj` objects always returns contiguous blocks, and `show_cp = TRUE` returns the first time index of each new block. Invalid `loss` values now raise an error.

## Other changes
- `posterior_estimate()` uses `maxZealousAttempts = 100` by default, avoiding spurious SALSO warnings when many change points are present.
- `posterior_estimate()` with `maxNClusters = 0` (default) sets no upper bound on the number of clusters, so SALSO no longer warns that the estimate "equals the default maximum possible number of clusters".
- `?posterior_estimate` now documents both methods and all available loss functions.
- GSL is located with `gsl-config` through a `configure` script, so the package installs from source on macOS and Linux with GSL installed in non-standard locations (e.g. Homebrew). Installation instructions for each operating system have been added to the README.
- Clarified the documentation of `prior_var_phi` (standard deviation of the proposal for `logit(phi)`) and `user_seed`.

# BayesChange 2.3.0

## New features
- Added argument `show_cp` to `posterior_estimate()` to provide also point estimation of the position of the change points in change point detection.

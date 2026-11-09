# BayesChange (development version)

## Bug fixes
- `user_seed` now fully determines the output of `detect_cp()` and `clust_cp()`. Previously only the GSL generator was seeded, while most draws came from R's random number generator, so two calls with the same `user_seed` could return different results. R's random number generator state is restored on exit; `user_seed = NULL` uses the current state of R's generator.
- Fixed the Metropolis-Hastings update of the autocorrelation parameter `phi` in `detect_cp()` for time series: the multivariate log-gamma function truncated its argument to an integer, the log full conditional evaluated to `-Inf` and every proposal was accepted. The same truncation also affected the integrated likelihood of multivariate time series.
- Fixed the integrated likelihood of univariate time series. It now equals the closed-form marginal of the model in Corradin, Danese & Ongaro (2022) (Eq. (11) with d = 1). Previously the function used by `clust_cp()` and in the split-merge steps had `log(1 - phi)` instead of `log(1 + phi)`, and blocks of length 1 and 2 were handled incorrectly.
- The update of `phi` in `detect_cp()` for univariate time series now uses the prior given by `a`, `b`, `c` (it previously used fixed hyperparameters), and its full conditional no longer omits a `log` in the `k_n` term (this also affected multivariate time series).
- Fixed the importance-sampling estimate of the normalisation constant used by `clust_cp()` (all kernels): the estimate was biased upwards by about `2 * log(B)` and used an incorrect proposal density, which made the algorithm favour too many clusters.
- Computational time is now recorded also when `print_progress = FALSE` (it was reported as 0).
- Random starting orders no longer contain empty blocks.
- `posterior_estimate()` for `DetectCpObj` objects always returns contiguous blocks, and `show_cp = TRUE` returns the first time index of each new block. Invalid `loss` values now raise an error.

## Other changes
- `posterior_estimate()` uses `maxZealousAttempts = 100` by default, avoiding spurious SALSO warnings when many change points are present.
- `?posterior_estimate` now documents both methods and all available loss functions.
- GSL is located with `gsl-config` through a `configure` script, so the package installs from source on macOS and Linux with GSL installed in non-standard locations (e.g. Homebrew). Installation instructions for each operating system have been added to the README.
- Clarified the documentation of `prior_var_phi` (standard deviation of the proposal for `logit(phi)`) and `user_seed`.

# BayesChange 2.3.0

## New features
- Added argument `show_cp` to `posterior_estimate()` to provide also point estimation of the position of the change points in change point detection.

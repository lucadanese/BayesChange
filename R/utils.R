# Seed R's RNG with `user_seed` for the duration of the calling function and
# restore the caller's RNG state on exit. Returns the seed to pass to C++.
#
# The C++ samplers draw from two generators: GSL (seeded directly with
# `user_seed`) and Armadillo (`arma::randu()`, `arma::randn()`, ...), which
# RcppArmadillo routes through R's RNG. Seeding both makes the output fully
# determined by `user_seed`, independently of any previous `set.seed()` call.
#
# If `user_seed` is NULL, R's RNG is left alone (so results follow
# `set.seed()`) and the GSL seed is drawn from it.
local_user_seed <- function(user_seed, envir = parent.frame()){
  if(is.null(user_seed)){
    return(sample.int(.Machine$integer.max, 1))
  }

  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if(had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL

  do.call(on.exit, list(substitute({
    if(had_seed){
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if(exists(".Random.seed", envir = globalenv(), inherits = FALSE)){
      rm(".Random.seed", envir = globalenv())
    }
  }, list(had_seed = had_seed, old_seed = old_seed)), add = TRUE), envir = envir)

  set.seed(user_seed)
  user_seed
}

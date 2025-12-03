# data-raw/epi_synthetic_multi.R
# Synthetic multivariate infection count dataset generated with
# the Doob–Gillespie algorithm

library(usethis)

# ---- 1. Prepare storage ----
data <- matrix(0, nrow = 3, ncol = 200)
inf_times <- vector("list", 3)

# ---- 2. Define beta vectors ----
betas <- list(
  c(rep(0.211, 120), rep(0.55, 80)),
  c(rep(0.215, 120), rep(0.52, 80)),
  c(rep(0.193, 30),  rep(0.53, 170))
)

# ---- 3. Simulate three epidemic processes ----
for(i in 1:3){
  # simulation of infection event times with sim_epi_data
  inf_times[[i]] <- sim_epi_data(
    S0 = 100000,
    I0 = 20,
    max_time = 200,
    beta_vec = betas[[i]],
    xi_0 = 1/8
  )

  # aggregate by integer time
  inf_times[[i]] <- table(floor(inf_times[[i]]))

  # fill row i of the data matrix
  data[i, as.numeric(names(inf_times[[i]]))] <- as.numeric(inf_times[[i]])
}

# ---- 4. Save dataset ----
epi_synthetic_multi <- data

usethis::use_data(epi_synthetic_multi, overwrite = TRUE)


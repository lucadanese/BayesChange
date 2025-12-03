# data-raw/epi_synthetic.R
# Script to generate the synthetic epidemiological dataset
# using a stochastic epidemic model (Doob–Gillespie algorithm)

library(usethis)
library(readr)

set.seed(3122025)

# ---- 1. Create synthetic beta vector ----
betas <- c(rep(0.2, 130), rep(0.55, 70))

# ---- 2. Simulate epidemic time series ----

# generate synthetic data with sim_epi_data()
inf_times <- sim_epi_data(
  S0 = 10000,
  I0 = 50,
  max_time = 200,
  beta_vec = betas,
  xi_0 = 1/8
)

# ---- 3. Aggregate infection events by integer time ----
inf_times_tab <- table(floor(inf_times))

inf_count <- matrix(0, nrow = 200, ncol = 1)
inf_count[as.numeric(names(inf_times_tab)), 1] <- as.numeric(inf_times_tab)

# ---- 4. Save as internal data object ----
epi_synthetic <- inf_count

usethis::use_data(epi_synthetic, overwrite = TRUE)


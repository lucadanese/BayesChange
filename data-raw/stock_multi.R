data_list_multi <- vector("list", length(companies))

for (i in seq_along(companies)) {

  df <- read_csv(
    paste0('data-raw/stock_data/', companies[i], ".csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      Mean = (High + Low) / 2,
      Date = as.Date(Date)
    ) %>%
    select(Date, Open, Mean, Close) %>%
    filter(Date >= First_Day, Date < Last_Day) %>%
    fill(Open, Close, Mean, .direction = "updown")

  # Standardize open/close
  open  <- scale(df$Open)[,1]
  close <- scale(df$Close)[,1]

  # Detrend each
  t <- seq_along(open)
  open  <- open  - fitted(lm(open  ~ t))
  close <- close - fitted(lm(close ~ t))

  # Combine into 2×T matrix
  M <- rbind(open, close)

  data_list_multi[[i]] <- M
}

# Convert to array 2 × 505 × 50
T_len <- ncol(data_list_multi[[1]])
stock_multi <- array(NA, dim = c(2, T_len, length(companies)))

for (i in seq_along(companies)) {
  stock_multi[, , i] <- data_list_multi[[i]]
}

dimnames(stock_multi) <- list(
  Price_Type = c("Open", "Close"),
  Time = paste0("t", seq_len(T_len)),
  Company = companies
)

# Save for package
usethis::use_data(stock_multi, overwrite = TRUE)

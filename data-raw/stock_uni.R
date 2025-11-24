library(dplyr)
library(tidyr)
library(readr)

First_Day <- as.Date("2020-01-01")
Last_Day  <- as.Date("2022-01-02")

companies <- c("MSFT","AAPL","AMZN","NVDA","GOOGL","GOOG","BRK-B","UNH","LLY",
               "JPM","XOM","AVGO","V","JNJ","PG","MA","HD","ADBE","COST","CVX",
               "MRK","PEP","KO","WMT","CRM","NFLX","ACN","MCD","BAC","LIN","AMD",
               "CSCO","TMO","INTC","ORCL","ABT","CMCSA","DIS","PFE","WFC","INTU",
               "VZ","PM","DHR","AMGN","IBM","QCOM","TXN","COP","UNP")

data_list <- vector("list", length(companies))

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

  # Extract Mean and transform to vector
  ts_vec <- df$Mean

  # Standardize
  ts_vec <- scale(ts_vec)[, 1]

  # Remove linear trend
  t <- seq_along(ts_vec)
  fit <- lm(ts_vec ~ t)
  ts_vec <- ts_vec - fitted(fit)

  data_list[[i]] <- ts_vec
}

# Convert to matrix (50 × 505)
stock_uni <- do.call(rbind, data_list)
rownames(stock_uni) <- companies

# Save for package
usethis::use_data(stock_uni, overwrite = TRUE)

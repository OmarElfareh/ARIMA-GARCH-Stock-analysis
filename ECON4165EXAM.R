# Import packages 
# Loading in packages 
library(forecast)
library(tidyverse)
library(tseries)
library(urca)
library(TSstudio)
library(readr)
library(ggplot2)
library(quantmod)
library(tseries)
library(urca)
library(xts)
library(fGarch)
library(rugarch)
library(readr)
library(xtable)
library(FinTS)

### DATA WORK AND ANALYSIS###


# Import data set
# Stock price data for Equinor 
getSymbols("EQNR", 
           src = "yahoo", 
           periodicity = "daily", 
           from = "2000/1/1",) 

# Making Equinor data into a time series 
EQNR_ts <- xts(EQNR$EQNR.Adjusted)

# Ploting Equinors time series 
plot(EQNR_ts, main = "EQNR Stock Price", ylab="Stock price in US $")

# Descriptive statistics 
summary(EQNR_ts)

# Plotting the log-adjusted time series
EQNR_log <- log(EQNR_ts)
plot(EQNR_log, main = "EQNR Stock Price (Log transformed)")
EQNR_log <- ts(EQNR_log)

# Checking if data is stationary and creating a LaTex table 
adf_result <- adf.test(EQNR_log)
adf_result

summary_table <- data.frame(
  TestStatistic = adf_result$statistic,
  PValue = adf_result$p.value,
  Lags = adf_result$parameter,
  row.names = "ADF Test"
)

# Convert the summary table to a LaTeX table
latex_table <- xtable(summary_table, caption = "ADF Test Results")

# Print the LaTeX table to the console
print(latex_table)

# difference the logged data
EQNR_returns <- diff(EQNR_log, lag = 1)
EQNR_returns <- na.omit(EQNR_returns)
plot(EQNR_returns, main = "EQNR stock returns", ylab = "Returns")

# Squaring the returns
EQNR_squared_returns <- EQNR_returns^2
plot(EQNR_squared_returns, main = "EQNR Squared Returns", ylab = "Squared Returns", xlab = "Time")


# Checking if data is stationary  
adf_result2 <- adf.test(EQNR_returns)
adf_result2
summary_table <- data.frame(
  TestStatistic = adf_result2$statistic,
  PValue = adf_result2$p.value,
  Lags = adf_result2$parameter,
  row.names = "ADF Test"
)

# Convert the summary table to a LaTeX table
latex_table <- xtable(summary_table, caption = "ADF Test Results 1. Difference")

# Print the LaTeX table to the console
print(latex_table)

# Perform the KPSS test on a time series
kpss.test(EQNR_returns)


# Differenced data correlogram
acf(EQNR_returns, main="Autocorrelation Function of EQNR returns")
pacf(EQNR_returns, main="Partial ACF of EQNR returns")


acf(EQNR_squared_returns, main="ACF of EQNR squared returns")
pacf(EQNR_squared_returns, main="PACF of EQNR squared returns")

### ARIMA MODEL ESTIMATION 


# Initialize an empty data frame to store AIC values and model orders
aic_table <- data.frame(
  p = integer(),
  d = integer(),
  q = integer(),
  AIC = numeric()
)

# Find optimal model AIC
best_model <- NULL
best_aic <- Inf

for (p in 0:5) {
  for (q in 0:5) {
    model <- arima(EQNR_returns, order = c(p, 0, q))
    aic <- AIC(model)
    if (aic < best_aic) {
      best_model <- model
      best_aic <- aic
    }
    # Add the model order and AIC to the table
    aic_table <- rbind(aic_table, data.frame(p = p, d = 0, q = q, AIC = aic))
  }
}

# Print the best model
print(best_model)

# Display the AIC table
print(aic_table)

# Create a LaTeX table for the AIC values
xtable(aic_table)

# Initialize an empty data frame to store BIC values and model orders
bic_table <- data.frame(
  p = integer(),
  d = integer(),
  q = integer(),
  BIC = numeric()
)

# Find optimal model BIC
best_model <- NULL
best_bic <- Inf

for (p in 0:5) {
  for (q in 0:5) {
    model <- arima(EQNR_returns, order = c(p, 0, q))
    bic <- BIC(model)
    if (bic < best_bic) {
      best_model <- model
      best_bic <- bic
    }
    # Add the model order and BIC to the table
    bic_table <- rbind(bic_table, data.frame(p = p, d = 0, q = q, BIC = bic))
  }
}

# Print the best model
print(best_model)

# Display the BIC table
print(bic_table)

# Create a LaTeX table for the BIC values
xtable(bic_table)


# Checking with auto.arima

ARIMA1 <- auto.arima(EQNR_returns, seasonal = TRUE)

# Extract model information
summary_ARIMA1 <- summary(ARIMA1)
summary_ARIMA1
model_order <- paste0("ARIMA(", paste(ARIMA1$order, collapse = ","), ")")
coef_summary <- summary_ARIMA1$coef
coef_se <- sqrt(diag(summary_ARIMA1$var.coef))

# Create a data frame with coefficients and standard errors
df_summary_ARIMA1 <- as.data.frame(cbind(coef_summary, coef_se))

# Add additional model information
model_info <- c(
  "Model" = model_order,
  "sigma^2" = summary_ARIMA1$sigma2,
  "Log Likelihood" = summary_ARIMA1$loglik,
  "AIC" = summary_ARIMA1$aic,
  "AICc" = summary_ARIMA1$aicc,
  "BIC" = summary_ARIMA1$bic
)

# Convert model_info to a data frame
df_model_info <- data.frame(
  coef_summary = model_info,
  coef_se = rep(NA, length(model_info)),
  row.names = names(model_info)
)

# Combine the coefficients and model_info data frames
df_complete1 <- rbind(df_summary_ARIMA1, df_model_info)

# Create LaTeX table
xtable(df_complete1)

# I also want the p-values 
# Extract coefficients and standard errors
arima_coefs <- coef(ARIMA1)
arima_se <- sqrt(diag(vcov(ARIMA1)))

# Calculate t-statistics and p-values
t_stats <- arima_coefs / arima_se
p_values <- 2 * (1 - pt(abs(t_stats), df = length(EQNR_returns) - length(arima_coefs)))

# Display results
results <- data.frame(
  Coefficient = names(arima_coefs),
  Estimate = arima_coefs,
  Std.Error = arima_se,
  TStatistic = t_stats,
  PValue = p_values
)

print(results)


# I also want to test ARIMA(1,0,1) due to the acf and pacf plot
# Extract model information
ARIMA2 <- Arima(EQNR_returns, order=c(1,0,1))
summary_ARIMA2 <- summary(ARIMA2)
model_order <- paste0("ARIMA(", paste(ARIMA2$order, collapse = ","), ")")
coef_summary <- summary_ARIMA2$coef
coef_se <- sqrt(diag(summary_ARIMA2$var.coef))

# Create a data frame with coefficients and standard errors
df_summary_ARIMA2 <- as.data.frame(cbind(coef_summary, coef_se))

# Add additional model information
model_info <- c(
  "Model" = model_order,
  "sigma^2" = summary_ARIMA2$sigma2,
  "Log Likelihood" = summary_ARIMA2$loglik,
  "AIC" = summary_ARIMA2$aic,
  "AICc" = summary_ARIMA2$aicc,
  "BIC" = summary_ARIMA2$bic
)

# Convert model_info to a data frame
df_model_info <- data.frame(
  coef_summary = model_info,
  coef_se = rep(NA, length(model_info)),
  row.names = names(model_info)
)

# Combine the coefficients and model_info data frames
df_complete2 <- rbind(df_summary_ARIMA2, df_model_info)

# Create LaTeX table
xtable(df_complete2)

# Pvalues for model 2 aswell
# Extract coefficients and standard errors
arima2_coefs <- coef(ARIMA2)
arima2_se <- sqrt(diag(vcov(ARIMA2)))

# Calculate t-statistics and p-values
t_stats2 <- arima2_coefs / arima2_se
p_values2 <- 2 * (1 - pt(abs(t_stats2), df = length(EQNR_returns) - length(arima2_coefs)))

# Display results
results2 <- data.frame(
  Coefficient = names(arima2_coefs),
  Estimate = arima2_coefs,
  Std.Error = arima2_se,
  TStatistic = t_stats2,
  PValue = p_values2
)

print(results2)


### DIAGNOSTICS CHECKS ###

Residuals_ARIMA1 <- residuals(ARIMA1, standardize = FALSE)

# Define the range of lag lengths to test
min_lag <- 5
max_lag <- 40
step <- 5

# Calculate the number of tests
n_tests <- (max_lag - min_lag) / step + 1

# Create a vector to store the p-values
p_values <- numeric(n_tests)

# Loop through the lag lengths and run the Ljung-Box test
counter <- 1
for (lag in seq(min_lag, max_lag, by = step)) {
  test_result <- Box.test(Residuals_ARIMA1, lag = lag, type = "Ljung-Box", fitdf=1)
  p_values[counter] <- test_result$p.value
  cat("Lag:", lag, "P-value:", test_result$p.value, "\n")
  counter <- counter + 1
}

# Print the p-values
p_values


Residuals_ARIMA2 <- residuals(ARIMA2, standardize = FALSE)

# Define the range of lag lengths to test
min_lag <- 5
max_lag <- 40
step <- 5

# Calculate the number of tests
n_tests <- (max_lag - min_lag) / step + 1

# Create a vector to store the p-values
p_values <- numeric(n_tests)

# Loop through the lag lengths and run the Ljung-Box test
counter <- 1
for (lag in seq(min_lag, max_lag, by = step)) {
  test_result <- Box.test(Residuals_ARIMA2, lag = lag, type = "Ljung-Box", fitdf=2)
  p_values[counter] <- test_result$p.value
  cat("Lag:", lag, "P-value:", test_result$p.value, "\n")
  counter <- counter + 1
}

# Print the p-values
p_values


# Histogram
hist(ARIMA1$residuals, main = "Histogram of Residuals for ARIMA(0,0,1)", 
     xlab = "Residuals", col = "lightblue", border = "black")

# Histogram
hist(ARIMA2$residuals, main = "Histogram of Residuals for ARIMA(1,0,1)", 
     xlab = "Residuals", col = "lightblue", border = "black")

# Q-Q plot
qqnorm(ARIMA1$residuals, main = "Q-Q Plot of Residuals for ARIMA(0,0,1)")
qqline(ARIMA1$residuals, col = "red")

# Q-Q plot
qqnorm(ARIMA2$residuals, main = "Q-Q Plot of Residuals for ARIMA(1,0,1)")
qqline(ARIMA2$residuals, col = "red")

# Check residuals 
checkresiduals(ARIMA1)
checkresiduals(ARIMA2)

acf(ARIMA1$residuals, main="ACF of ARIMA(0,0,1) residuals")
pacf(ARIMA1$residuals, main="PACF of ARIMA(0,0,1)residuals")

acf(ARIMA2$residuals, main="ACF of ARIMA(1,0,1) residuals")
pacf(ARIMA2$residuals, main="PACF of ARIMA(1,0,1)residuals")

# Out-of-sample forecasting

# Splitting the time series into training and test set
split_EQNR <- ts_split(EQNR_returns, sample.out = 500)
training <- split_EQNR$train
testing <- split_EQNR$test

# Out-of-sample forecast for ARIMA1
forecast <- vector()
for (i in 1:length(testing)){
  r.temp <- EQNR_returns[(1:(length(training)-1+i))]
  arima.temp <- Arima(r.temp, order = c(0,0,1))
  forecast[i] <- predict(arima.temp,1)$pred
}

# Out-of-sample forecast for ARIMA2
forecast2 <- vector()
for (i in 1:length(testing)){
  r.temp2 <- EQNR_returns[(1:(length(training)-1+i))]
  arima.temp2 <- Arima(r.temp2, order = c(1,0,1))
  forecast2[i] <- predict(arima.temp2,1)$pred
}

# Plot actual returns and forecasts from both models
plot(as.numeric(testing), type = "l", col = "black", 
     main = "Actual Returns & Forecasts ARIMA(0,0,1) vs ARIMA(1,0,1)", ylab = "Returns")
lines(forecast, col = "red")
lines(forecast2, col = "blue")
legend("topright", legend = c("Actual Returns", "Forecast ARIMA(0,0,1)", "Forecast ARIMA(1,0,1)"), 
       col = c("black", "red", "blue"), lty = 1, cex = 0.8)

# Error Metrics 

# RMSE and MAE for ARIMA1 forecasts 
RMSE <- sqrt(mean((testing - forecast)^2))
MAE <- mean(abs(testing - forecast))
cat("RMSE:", RMSE, "\nMAE:", MAE)

# RMSE and MAE for ARIMA2 forecasts 
RMSE <- sqrt(mean((testing - forecast2)^2))
MAE <- mean(abs(testing - forecast2))
cat("RMSE:", RMSE, "\nMAE:", MAE)

# MAPE for ARIMA 1 forecasts 
MAPE <- mean(abs((testing - forecast) / testing)) * 100
cat("MAPE:", MAPE, "%")

# MAPE for ARIMA 2 forecasts 
MAPE <- mean(abs((testing - forecast2) / testing)) * 100
cat("MAPE:", MAPE, "%")

#sMAPE for ARIMA 1 forecasts 
sMAPE <- mean(2 * abs(testing - forecast) / (abs(testing) + abs(forecast))) * 100
cat("sMAPE:", sMAPE, "%")

#sMAPE for ARIMA 2 forecasts 
sMAPE <- mean(2 * abs(testing - forecast2) / (abs(testing) + abs(forecast2))) * 100
cat("sMAPE:", sMAPE, "%")

# Arch Test for ARIMA1
ArchTest(ARIMA1$residuals)

# Arch Test for ARIMA2
ArchTest(ARIMA2$residuals)

"Given the very low p-value, we can reject the null hypothesis that there are 
no ARCH effects in the EQNR_returns data. This suggests that the data exhibits
time-varying volatility, and a model like GARCH might be more appropriate to 
capture the volatility structure in the data."

### GARCH MODEL ESTIMATION

# Find the best ARMA(0,1)-GARCH(p,q) model using a loop with the BIC information criteria

Result_BIC <- data.frame(Model="m",BIC=0)

for (i in 1:4){
  for (j in 1:4){
    fit=garchFit(substitute(~arma(0,1)+garch(p,q),list(p=i,q=j)),data = EQNR_returns,trace=F)
    Result_BIC=rbind(Result_BIC,data.frame(Model=paste("m-",i,"-",j),BIC=fit@fit$ics[2]))
  }
}

view(Result_BIC)
Result_BIC_xt <- xtable(Result_BIC)
print.xtable(Result_BIC_xt, caption.placement = "top", caption = "Best ARMA(0,1)-GARCH(p,q) Models based on BIC", label = "table:garch_bic", table.placement = "ht", include.rownames = FALSE, digits = 10)


garch1 <- garchFit(data = EQNR_returns, formula = ~ arma(0,1) + garch(1,1), trace = F) 
garch1

### Diagnostics checks 

# Ljung Box Test on the standardized residuals of the ARMA(0,1)-GARCH(1,1) model
std_residualsgarch <- residuals(garch1, standardize = TRUE)

# Define the range of lag lengths to test
min_lag <- 5
max_lag <- 40
step <- 5

# Calculate the number of tests
n_tests <- (max_lag - min_lag) / step + 1

# Create a vector to store the p-values
p_values <- numeric(n_tests)

# Loop through the lag lengths and run the Ljung-Box test
counter <- 1
for (lag in seq(min_lag, max_lag, by = step)) {
  test_result <- Box.test(std_residualsgarch, lag = lag, type = "Ljung-Box", fitdf = 3)
  p_values[counter] <- test_result$p.value
  cat("Lag:", lag, "P-value:", test_result$p.value, "\n")
  counter <- counter + 1
}

# Print the p-values
p_values


# Ljung Box Test on the squared standardized residuals of the ARMA(0,1)-GARCH(1,1) model
squared_std_residuals <- std_residualsgarch^2

# Define the range of lag lengths to test
min_lag <- 5
max_lag <- 40
step <- 5

# Calculate the number of tests
n_tests <- (max_lag - min_lag) / step + 1

# Create a vector to store the p-values
p_values <- numeric(n_tests)

# Loop through the lag lengths and run the Ljung-Box test
counter <- 1
for (lag in seq(min_lag, max_lag, by = step)) {
  test_result <- Box.test(squared_std_residuals, lag = lag, type = "Ljung-Box", fitdf = 3)
  p_values[counter] <- test_result$p.value
  cat("Lag:", lag, "P-value:", test_result$p.value, "\n")
  counter <- counter + 1
}

# Print the p-values
p_values

# Out-of-Sample forecast
predictions <- vector(mode="numeric", length=0)
res <- residuals(garch1)

for (i in 1:500){
  
  r.temp <- EQNR_returns[1:(length(EQNR_returns)-500+i)]  
  garch2 <- garchFit(data = r.temp, formula = ~ arma(0,0) +garch(1,1), trace = F) 
  
  predictions[i] <- predict(garch2, 1)$standardDeviation[1]^2
}
length(predictions)
plot.ts(res[(length(EQNR_returns)-500):length(EQNR_returns)]^2, 
        main = "Actual Squared Residuals & Forecasts ARMA(0,0)-GARCH(1,1)",
        ylab = "Squared Residuals",
        xlab = "Time")
lines(predictions, col = "red")

# Adding a legend to the plot
legend("topright", 
       legend = c("Actual Squared Residuals", "Squared Residuals Forecasts"), 
       col = c("black", "red"), 
       lty = 1, 
       cex = 0.8)

# Forecast error metrics 
# Actual squared residuals for the last 500 data points
actual_squared_residuals <- res[(length(EQNR_returns)-500):length(EQNR_returns)]^2

# Calculate RMSE
RMSE <- sqrt(mean((actual_squared_residuals - predictions)^2))
cat("RMSE:", RMSE, "\n")

# Calculate MAE
MAE <- mean(abs(actual_squared_residuals - predictions))
cat("MAE:", MAE, "\n")

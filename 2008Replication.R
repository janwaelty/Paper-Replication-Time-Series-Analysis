###############################################################################
#   Replication Paper 2008: Are combination forecasts of S&P 500 
#                           volatility statistically superior?
#   Jan Wälty 20-704-227
###############################################################################
library(tidyverse)
library(stargazer)
library(tsbox)
library(rugarch)
library(forecast)
library(tseries)
library(yfR)
library(fGarch)
library(reshape2)
library(dplyr)

###############################################################################
#data pre-processing
###############################################################################

#fetch  S&P 500 as well as VIX data using yfR package 
first_date <-as.Date("1989-12-29")
last_date <- as.Date("2003-10-18")
freq <- "daily"
data_SP <- yf_get(tickers = '^GSPC' , first_date = first_date, 
                  last_date = last_date, freq = "daily")
data_VIX <- yf_get(tickers = "^VIX", first_date = first_date, 
                   last_date = last_date, freq = "daily")

# exclude first row 
data_SP <- data_SP[-1, ]

# choose relevant columns
data_SP <- data_SP %>%
  select(ref_date,ret_closing_prices, price_low, price_high)
data_VIX <- data_VIX %>%
  select(ref_date, price_close)

# plot return series
plot_data<- data_SP %>% select(ref_date, ret_closing_prices) 
png(file = "C:/Universität/EWF/return_SP.png", height = 7, width = 7, 
    units = "in", res = 800)
ggplot(plot_data, aes(x = ref_date, y = ret_closing_prices)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Return"
  )
dev.off()

# de-mean return function
demean = function(data){
  return(data-mean(data))
}
# function for effective volatility (for 22 days average) using squared returns
# use square root of mean of demeaned squared returns
realized_vol_y_squared = function(start){
  return(sqrt(mean(demean(data_SP$ret_closing_prices[start:(start + 21)])^2)))
}

#function for effective volatility (for 22 days average) using high/low ratio
realized_vol_high_low = function(start){
  return(mean((data_SP$price_high[start:(start + 21)] - 
                 data_SP$price_low[start:(start + 21)])/(data_SP$price_low[start:(start + 21)])))
}

###############################################################################
#implement functions for different models
###############################################################################
#Notice : For Garch family models, returns do not need to be demeaned as model 
#         implementation internally does that

#i) standard VIX model

vix <- function(start){
  return(data_VIX$price_close[start] / sqrt(252) * 0.01)
  'VIX yields annualized expected standard deviation in percent, so 
  have to convert to daily standard deviation in decimal notation 
  (252 trading days)'
}

#ii) GARCH(1,1)

garch <- function(start, end){
  arch.fit <- garchFit(~garch(1,1), data = data_SP$ret_closing_prices[start:end])
  forecast <- predict(arch.fit, n.ahead = 22)
  return(mean(forecast$standardDeviation))
  # returns mean standard deviation forecast over 22 days
}

#iii) GJR-Garch(1,1)

gjr_garch <- function(start, end){
  # rugarch packages needed; hybrid solver required for convergence
  spec= ugarchspec(variance.model = list(model="gjrGARCH",garchOrder = c(1, 1))) 
  garch_fit <- ugarchfit(spec = spec,
                         data = data_SP$ret_closing_prices[start:end],
                         solver = 'hybrid')
  garch_forecast <- ugarchforecast(garch_fit, n.ahead = 22)
  return(mean(garch_forecast@forecast$sigmaFor)) 
  # returns mean standard deviation forecast over 22 days
}

#iv) ARCH(1)

arch <- function(start, end){
  arch.fit <- garchFit(~garch(1,0), data = data_SP$ret_closing_prices[start:end])
  forecast <- predict(arch.fit, n.ahead = 22)
  return(mean(forecast$standardDeviation))
}

###############################################################################
# rolling window procedure
###############################################################################

#number of prediction intervals (1000 observations each)
num_predictions = 2460

# lists for volatility predictions
vol_pred_vix = list()
vol_pred_garch = list()
vol_pred_gjr = list()
vol_pred_arch = list()

# get volatility estimates from models 
for(idx in 1:num_predictions){
  vol_pred_vix[idx] <- vix((idx+ 999))
  vol_pred_garch[idx] <- garch(idx, (idx+999))
  vol_pred_gjr[idx] <- gjr_garch(idx, (idx+999))
  vol_pred_arch[idx] <- arch(idx, (idx+999))
}
###############################################################################
# Squared return proxy calculations 

MSE_vix_list = list()
MSE_garch_list = list()
MSE_gjr_list = list()
MSE_arch_list = list()
squared_list = list()


for(idx in 1: num_predictions){
  #proxy for specific horizon
  squared_list[idx] <- realized_vol_y_squared((idx+1000))
  #squared deviation to high_low proxy
  MSE_vix_list[idx] <- (vol_pred_vix[[idx]]- 
                          realized_vol_y_squared((idx+1000)))^2
  MSE_garch_list[idx] <- (vol_pred_garch[[idx]]-
                            realized_vol_y_squared((idx+1000)))^2
  MSE_gjr_list[idx] <- (vol_pred_gjr[[idx]]- 
                          realized_vol_y_squared((idx+1000)))^2
  MSE_arch_list[idx] <- (vol_pred_arch[[idx]]- 
                           realized_vol_y_squared((idx+1000)))^2
}

# sum over all deviations (without averaging as in paper)
MSE_vix <- sum(unlist(MSE_vix_list))
MSE_garch <- sum(unlist(MSE_garch_list))
MSE_gjr <- sum(unlist(MSE_gjr_list))
MSE_arch <- sum(unlist(MSE_arch_list))


###############################################################################
# high-low ratio calculations 

MSE_vix_list_high_low = list()
MSE_garch_list_high_low  = list()
MSE_gjr_list_high_low  = list()
MSE_arch_list_high_low  = list()
squared_list_high_low = list()

for(idx in 1: num_predictions){
  squared_list_high_low[idx] <- realized_vol_high_low((idx+1000))
  MSE_vix_list_high_low[idx] <- (vol_pred_vix[[idx]]- 
                                   realized_vol_high_low((idx+1000)))^2
  MSE_garch_list_high_low[idx] <- (vol_pred_garch[[idx]]-
                                     realized_vol_high_low((idx+1000)))^2
  MSE_gjr_list_high_low[idx] <- (vol_pred_gjr[[idx]]- 
                                   realized_vol_high_low((idx+1000)))^2
  MSE_arch_list_high_low[idx] <- (vol_pred_arch[[idx]]- 
                                    realized_vol_high_low((idx+1000)))^2
}

MSE_vix_high_low <- sum(unlist(MSE_vix_list_high_low))
MSE_garch_high_low <- sum(unlist(MSE_garch_list_high_low))
MSE_gjr_high_low <- sum(unlist(MSE_gjr_list_high_low))
MSE_arch_high_low <- sum(unlist(MSE_arch_list_high_low))

###############################################################################
# Combination forecast (using squared returns as proxy)

#prepare data
data_to_predict <- data.frame(squared_list = unlist(squared_list),
                              vol_pred_vix = unlist(vol_pred_vix), 
                              vol_pred_garch = unlist(vol_pred_garch),
                              vol_pred_gjr = unlist(vol_pred_gjr),
                              vol_pred_arch = unlist(vol_pred_arch))

# combine all volatility models and find weights by multiple regression 
# Notice: Slightly differs from implementation in paper that re-estimates the 
#         weights (but requires intraday data)
combined_model <- lm(data_to_predict$squared_list ~ data_to_predict$vol_pred_vix + 
                       + data_to_predict$vol_pred_garch + 
                       data_to_predict$vol_pred_gjr + 
                       data_to_predict$vol_pred_arch)

#predict forecasts based on weights from above
combined_forecast <- predict(combined_model, newdata = data_to_predict)
MSE_combined_forecast_squared <- list()
MSE_combined_forecast_high_low <- list()

for(idx in 1: num_predictions){
  # compute deviations directly for both proxies
  MSE_combined_forecast_squared[idx] <- (combined_forecast[[idx]]- 
                                           realized_vol_y_squared((idx+1000)))^2
  MSE_combined_forecast_high_low[idx] <- (combined_forecast[[idx]]- 
                                            realized_vol_high_low((idx+1000)))^2
}

#sum over all deviations (without averaging)
MSE_combined_squared <- sum(unlist(MSE_combined_forecast_squared))
MSE_combined_high_low <- sum(unlist(MSE_combined_forecast_high_low))

###############################################################################
#Summary of results

#print MSE table (with squared returns as proxy for realized volatility)
table_squared <- data.frame(
  "Garch" = round(MSE_garch,3),
  "GJR-Garch" = round(MSE_gjr,3),
  "ARCH" = round(MSE_arch,3),
  "VIX" = round(MSE_vix,3), 
  "Combined" = round(MSE_combined_squared,3)
)
print(table_squared)


#print MSE table (with high-low ratio as proxy for realized volatility)
table_high_low <- data.frame(
  "Garch" = round(MSE_garch_high_low,3),
  "GJR-Garch" = round(MSE_gjr_high_low,3),
  "ARCH" = round(MSE_arch_high_low,3),
  "VIX" = round(MSE_vix_high_low,3),
  "Combined" = round(MSE_combined_high_low, 3)
)
print(table_high_low)

###############################################################################
# plot deviation series against squared returns series
###############################################################################

# transform data to vectors
vol_pred_vix <- unlist(vol_pred_vix)
vol_pred_garch <- unlist(vol_pred_garch)
vol_pred_gjr <- unlist(vol_pred_gjr)
vol_pred_arch <- unlist(vol_pred_arch)
squared_list <- unlist(squared_list)
squared_list_high_low <- unlist(squared_list_high_low)


png(file = "C:/Universität/EWF/multiple.png", height = 7, width = 7, 
    units = "in", res = 1000)

par(mfrow = c(2, 3), oma = c(2, 2, 0, 0))   

# Plot 1: Squared returns against Arch forecast
plot(squared_list, type = "l", col = "black", lwd = 2,
     main = "RV and ARCH",
     xlab = "Prediction horizon",
     ylab = "Standard deviation")
lines(vol_pred_arch, col = "red", lwd = 1.5)

# Plot 2: Squared returns against and GJR-Garch forecast
plot(squared_list, type = "l", col = "black", lwd = 2,
     main = "RV and GJR",
     xlab = "Prediction horizon",
     ylab = "Standard deviation")
lines(vol_pred_gjr, col = "blue", lwd = 1.0)

# Plot 3: Squared returns against Garch forecast
plot(squared_list, type = "l", col = "black", lwd = 2,
     main = "RV and GARCH",
     xlab = "Prediction horizon",
     ylab = "Standard deviation")
lines(vol_pred_garch, col = "green", lwd = 1.0)

# Plot 4: Squared returns against VIX forecast
plot(squared_list, type = "l", col = "black", lwd = 2,
     main = "RV and VIX",
     xlab = "Prediction horizon",
     ylab = "Standard deviation")
lines(vol_pred_vix, col = "orange", lwd = 1.5)

# Plot 5: Squared returns against combination forecast
plot(squared_list, type = "l", col = "black", lwd = 2, 
     main = "RV and combined forecast",
     xlab = "Prediction horizon",
     ylab = "Standard deviation")
lines(combined_forecast, col = "violet", lwd = 1.0)

dev.off()


###############################################################################
# Plot RV models  against time series of returns 
###############################################################################

returns_pred_time <- data_SP$ret_closing_prices[1001:3460]
png(file = "C:/Universität/EWF/rv.png", height = 7, width = 7,
    units = "in", res = 800)
ylim <- c(-0.06, 0.08)
par(mfrow = c(2, 1))
plot(returns_pred_time, type = 'l', col = "blue", lwd = 2,
     main = "RV (squared returns) and return series",
     xlab = "Prediction horizon",
     ylab = "Return / SD", ylim = ylim)
lines(squared_list, col = "black", lwd = 3)
legend("topright", legend = c("RV", "S&P 500 returns"), 
       col = c("black", "blue"), lty = 1, lwd = c(2, 2), cex = 0.8)

plot(returns_pred_time, type = 'l', col = "blue", lwd = 2,
     main = "RV (high-low deviation) and returns",
     xlab = "Observation number",
     ylab = "Return / SD", ylim = ylim)
lines(squared_list_high_low, col = "black", lwd = 3)
legend("topright", legend = c("RV", "S&P 500 returns"), 
       col = c("black", "blue"), lty = 1, lwd = c(2, 2), cex = 0.8)

dev.off()





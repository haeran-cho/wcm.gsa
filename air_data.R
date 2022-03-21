## Analysis of nitrogen oxides concentrations in London

source('main.R')

no <- read.csv("AirQualityData.csv", header = TRUE) 
# contains dates and daily NO2 and NOx concentrations
# measured at Marylebone Road in London, UK,
# from September 1, 2000 to September 30, 2020

hols <- read.csv("ukbankholidays.csv", header = TRUE)
# contains all UK bank holidays

allhols <- as.Date(hols$UK.BANK.HOLIDAYS, "%d-%b-%Y")
day <- rep(1:7, ceiling(length(no$NO2)/7))[1:length(no$NO2)]
no$day <- day
no <- no[!is.na(no$NO2), ]
dates <- as.Date(no$Date)

## for analysing NO2 concentrations
xx <- sqrt(no$NO2)

## for analysing NOX concentrations
xx <- sqrt(no$NOX)

## seasonality and trend removal

dat <- data.frame(months = format.Date(dates, "%m"), days = factor(no$day), 
                  is.hol = dates %in% allhols, y = xx)

fit <- lm(y ~ ., data = dat[format.Date(dates, "%Y") < "2011" &
                                    format.Date(dates, "%Y") > "2003", ])

x <- xx - predict(fit, newdata = dat)

## change point analysis
w <- wcm.gsa(x, max.iter = 10, min.len = 10)
dates[w$cp]

plot(x, type = 'l', xaxt = 'n', xlab = '', ylab = '')
axis(1, at = seq(1, length(x), length.out = 20), label = format.Date(dates, "%Y-%m")[seq(1, length(x), length.out = 20)])
abline(v = w$cp, col = 2)

fhat <- x * 0
brks <- c(0, w$cp, length(x))
for(ii in 1:(length(w$cp) + 1)){
  int <- (brks[ii] + 1):brks[ii + 1]
  fhat[int] <- mean(x[int])
}

acf(x)
acf(x - fhat)


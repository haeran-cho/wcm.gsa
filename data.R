## Air Quality data analysis

no <- read.csv("AirQualityData.csv", header = TRUE) 
hols <- read.csv("ukbankholidays.csv", header = TRUE)

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
w <- wem.gsc(x, max.iter = 10, double.cusum = !TRUE)
dates[w$cp]
dates[w$rcp]

plot(x, type = 'l', xaxt = 'n', xlab = '', ylab = '')
axis(1, at = seq(1, length(x), length.out = 20), label = format.Date(dates, "%Y-%m")[seq(1, length(x), length.out = 20)])
abline(v = w$rcp, col = 2)

brks <- c(0, w$rcp, length(x)); fhat <- x*0
for(kk in 1:(length(brks) - 1)) fhat[(brks[kk] + 1):brks[kk + 1]] <- mean(x[(brks[kk] + 1):brks[kk + 1]])
lines(fhat, col = 4, lwd = 2)

acf(x)
acf(x - fhat)

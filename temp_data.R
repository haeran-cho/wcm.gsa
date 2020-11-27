## Analysis of Hadley Centre central England temperature (HadCET) data 

source('main.R')
dat <- read.csv('hadcet.csv')
# contains 4 columns including years and 
# yearly average of the monthly mean, maximum and minimum temperatures from 1878 to 2019


x <- dat$mean
# x <- dat$max
# x <- dat$min

## change point analysis
w <- wem.gsc(x, min.len = 10, p.max = 5, double.cusum = !TRUE)
dat$dates[w$cp]
dat$dates[w$rcp]

plot(x, type = 'l', xaxt = 'n', xlab = '', ylab = '')
axis(1, at = seq(1, length(x), length.out = 20), label = dat$dates[seq(1, length(x), length.out = 20)])
abline(v = w$rcp, col = 2)

brks <- c(0, w$rcp, length(x)); fhat <- x*0
for(kk in 1:(length(brks) - 1)) fhat[(brks[kk] + 1):brks[kk + 1]] <- mean(x[(brks[kk] + 1):brks[kk + 1]])
lines(fhat, col = 4, lwd = 2)

acf(x)
acf(x - fhat)

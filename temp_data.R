## Analysis of Hadley Centre central England temperature (HadCET) data 

source('main.R')
dat <- read.csv('hadcet.csv')
# contains 4 columns including years and 
# yearly average of the monthly mean, maximum and minimum temperatures from 1878 to 2019

x <- dat$mean
x <- dat$max
x <- dat$min

matplot(dat[, 2:4], type = 'l', xaxt = 'n', xlab = '', ylab = '')
axis(1, at = seq(1, length(x), length.out = 20), label = dat$dates[seq(1, length(x), length.out = 20)])


## change point analysis
w <- wcm.gsa(x, min.len = 10, p.max = 5)
dat$dates[w$cp]

plot(x, type = 'l', xaxt = 'n', xlab = '', ylab = '')
axis(1, at = seq(1, length(x), length.out = 20), label = dat$dates[seq(1, length(x), length.out = 20)])
abline(v = w$cp, col = 2)

fhat <- x * 0
brks <- c(0, w$cp, length(x))
for(ii in 1:(length(w$cp) + 1)){
  int <- (brks[ii] + 1):brks[ii + 1]
  fhat[int] <- mean(x[int])
}

acf(x)
acf(x - fhat)

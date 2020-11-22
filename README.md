# wem.gsc
Software accompanying 
> H. Cho and P. Fryzlewicz (2020) "Multiple change point detection under serial dependence: wild energy maximisation and gappy Schwarz criterion".

The main routine of WEM.gSC is contained in main.R (`wem.gsc`) and the code used for air quality data analysis is contained in data.R.

To use WEM.gSC, do the following:

- Source main.R into `R`.
- Read the descriptions within main.R.

For example,

```{r}
source("main.R")
set.seed(111)
f <- rep(c(0, 5, 2, 8, 1, -2), c(100, 200, 200, 50, 200, 250))
x <- f + arima.sim(list(ar = c(.75, -.5), ma = c(.8, .7, .6, .5, .4, .3)), n = length(f), sd = 1)
wem.gsc(x, double.cusum = TRUE)
```

If you have any questions, contact haeran.cho@bristol.ac.uk

# wcm.gsa
Software accompanying 
> H. Cho and P. Fryzlewicz (2021) ["Multiple change point detection under serial dependence: wild contrast maximisation and gappy Schwarz algorithm".](https://arxiv.org/abs/2011.13884)

- The main routine of WCM.gSa (`wcm.gsa`) is contained in main.R. 
- Codes used for analysing nitrogen oxides concentrations in London (contained in `AirQualityData.csv`) are provided in air_data.R.
- Codes used for analysing Hadley Centre central England temperature data (contained in `hadcet.csv`) are provided in temp_data.R.

To use WCM.gSa, do the following:

- Source main.R into `R`.
- Read the descriptions for `wcm.gsa` within main.R.

For example,

```{r}
source("main.R")
set.seed(111)
f <- rep(c(0, 5, 2, 8, 1, -2), c(100, 200, 200, 50, 200, 250))
x <- f + arima.sim(list(ar = c(.75, -.5), ma = c(.8, .7, .6, .5, .4, .3)), n = length(f), sd = 1)
wcm.gsa(x)
```

If you have any questions, contact haeran.cho@bristol.ac.uk

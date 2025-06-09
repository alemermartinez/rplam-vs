In this package, a robust estimation and variable selection procedure for partially linear additive models is performed using the real dataset from nutritional epidemiology.

Let's first load the script.

``` r
source("R/rplam-vs-fn.R")
```

Let's begin by reading the data.

``` r
datos <- read.table("Plasma_Retinol.txt", sep="\t")
str(datos)
```

Since a robust approach will be applied, we standarizing the continuos covariates robustly.

``` r
age <- (datos$V1 - median(datos$V1))/mad(datos$V1)
sex <- datos$V2-1
smokstat <- datos$V3
smok1 <- as.numeric(datos$V3==2)
smok2 <- as.numeric(datos$V3==3)
quetelet <- (datos$V4 - median(datos$V4))/mad(datos$V4)
vituse <- datos$V5
vit1 <- as.numeric(datos$V5==1)
vit2 <- as.numeric(datos$V5==2)
calories <- (datos$V6 - median(datos$V6))/mad(datos$V6)
fat <- (datos$V7-median(datos$V7))/mad(datos$V7)
fiber <- (datos$V8-median(datos$V8))/mad(datos$V8)
alcohol <- (datos$V9-median(datos$V9))/mad(datos$V9)
cholesterol <- (datos$V10-median(datos$V10))/mad(datos$V10)
betadiet <- (datos$V11-median(datos$V11))/mad(datos$V11)
retdiet <- (datos$V12-median(datos$V12))/mad(datos$V12)
retplasma <- (datos$V14-median(datos$V14))/mad(datos$V14)
Z <- cbind(sex,smok1,smok2,quetelet,vit1,vit2,calories,
           fat,alcohol,betadiet)
X <- cbind(age,cholesterol,fiber)
```

And also the response variable.

``` r
betaplasma <- (datos$V13-median(datos$V13))/mad(datos$V13)
y <- betaplasma
```

It can be appreciated an extreme outlier in 'alcohol'.
``` r
boxplot(alcohol)
```
![](README_files/figure-markdown_github/alcohol-1.png)

We will apply a robust estimator that simultaneoulsy select variables from the linear 
and the additive parts of the model.

We consider two grids for the auxiliary parameters.
``` r
grid.la1 <- seq(0, 0.1, by=0.01)
grid.la2 <- seq(0, 2, by=0.1)
```
and compute the proposal using cubic splines. Take into account that the computation
of the estimator over these grid takes about 723 secs in an Intel Core i7-10700 CPU @ 2.90GHz × 16.
``` r
degree.spline <- 3
system.time(
fit.rob <- plam.rob.vs(y, Z, X, degree.spline=degree.spline, grid.la1=grid.la1, grid.la2=grid.la2)
)
```

The results obtained are the following:
```
fit.rob$nknots
fit.rob$coef.const
fit.rob$coef.lin
fit.rob$la1
fit.rob$la2
```

Are there covariables irrelevant for the model?
``` r
fit.rob$is.zero
```
It can be appreciated that only five covariates are considered important for the model for the 
robust approach.

Let's see if the proposal identifies any large residuals:
``` r
res <- y-fit.rob$prediction
aa <- boxplot(y-fit.rob$prediction, col="lightblue")
length(aa$out)
```
17 observations were identified by the boxplot. These observations corresponds to observations
``` r
in.ro <- (1:length(res))[ res %in% aa$out ]
in.ro
```

In this package, a robust estimation and variable selection procedure
for partially linear additive models is performed using the real dataset
from nutritional epidemiology.

Let’s first load the script.

``` r
source("R/rplam-vs-fn.R")
```

Let’s begin by reading the data.

``` r
datos <- read.table("Plasma_Retinol.txt", sep="\t")
str(datos)
```

    ## 'data.frame':    315 obs. of  14 variables:
    ##  $ V1 : int  64 76 38 40 72 40 65 58 35 55 ...
    ##  $ V2 : int  2 2 2 2 2 2 2 2 2 2 ...
    ##  $ V3 : int  2 1 2 2 1 2 1 1 1 2 ...
    ##  $ V4 : num  21.5 23.9 20 25.1 21 ...
    ##  $ V5 : int  1 1 2 3 1 3 2 1 3 3 ...
    ##  $ V6 : num  1299 1032 2372 2450 1952 ...
    ##  $ V7 : num  57 50.1 83.6 97.5 82.6 56 52 63.4 57.8 39.6 ...
    ##  $ V8 : num  6.3 15.8 19.1 26.5 16.2 9.6 28.7 10.9 20.3 15.5 ...
    ##  $ V9 : num  0 0 14.1 0.5 0 1.3 0 0 0.6 0 ...
    ##  $ V10: num  170.3 75.8 257.9 332.6 170.8 ...
    ##  $ V11: int  1945 2653 6321 1061 2863 1729 5371 823 2895 3307 ...
    ##  $ V12: int  890 451 660 864 1209 1439 802 2571 944 493 ...
    ##  $ V13: int  200 124 328 153 92 148 258 64 218 81 ...
    ##  $ V14: int  915 727 721 615 799 654 834 825 517 562 ...

Since a robust approach will be applied, we standarizing the continuos
covariates robustly.

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

It can be appreciated an extreme outlier in ‘alcohol’.
![](README_files/figure-markdown_github/alcohol-1.png)

We will apply a robust estimator that simultaneoulsy select variables
from the linear and the additive parts of the model.

We consider two grids for the auxiliary parameters.

``` r
grid.la1 <- seq(0, 0.1, by=0.01)
grid.la2 <- seq(0, 2, by=0.1)
```

and compute the proposal using cubic splines. Take into account that the
computation of the estimator over these grid takes about 723 secs in an
Intel Core i7-10700 CPU @ 2.90GHz × 16.

``` r
degree.spline <- 3
system.time(
fit.rob <- plam.rob.vs(y, Z, X, degree.spline=degree.spline, grid.la1=grid.la1, grid.la2=grid.la2)
)
```

    FALSE [1] 0
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30351e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.22094e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.7241e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72431e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.15322e-20
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.19705e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.10562e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.93798e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.87856e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.3975e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.38424e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.10256e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32156e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.12366e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63906e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.57562e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08093e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.12606e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.78327e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.03877e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.71752e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.08031e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20633e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.23101e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.48303e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.29789e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.88428e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.98374e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.31733e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.71458e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.40466e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34624e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.00323e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.2897e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.04547e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.00359e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.08431e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.84387e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.4317e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.72282e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.40925e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.81341e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.44652e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.23604e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27087e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.54394e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.6237e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4034e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.81052e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.63088e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00451e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01083e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.70571e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.55277e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.18231e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.07864e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.3054e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.58661e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26029e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01339e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.99757e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65633e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.90634e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98159e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.8351e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.75673e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.31891e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36023e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18041e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34999e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.57477e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77694e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76272e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.39501e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4437e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.75251e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.12519e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.37075e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.79076e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.13536e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.15004e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74618e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.88264e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.8799e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08502e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.66711e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4433e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.37642e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.10934e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18906e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.1518e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63295e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36521e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.52928e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25396e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.57438e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.44866e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.67141e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56171e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.44382e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58843e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.80632e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.25512e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03463e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01293e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.60189e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.52225e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.64483e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58426e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.31275e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.6041e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47818e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19154e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08057e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28012e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.26271e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49788e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.80991e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08172e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.97768e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.00675e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.45511e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32367e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13112e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.70335e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16687e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14854e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.97849e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.64589e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.4679e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02364e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.81515e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76389e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.70002e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33852e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.84172e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67812e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.98459e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16471e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.61915e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65046e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.7743e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58521e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.25634e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36996e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.06475e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59273e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.28936e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.67829e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.33034e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.89594e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72851e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.68449e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.85906e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4003e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.12469e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.23712e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.43685e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73344e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74001e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.75617e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34548e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.97069e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.93291e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14857e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.55878e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25775e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.25303e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.35393e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98151e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.94617e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63507e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11891e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.43016e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.153e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.24299e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00226e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.70463e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.78318e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.15204e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.40108e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59033e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08696e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.06498e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.1886e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.28219e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01571e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.82495e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30959e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.7131e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.82519e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.23954e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.27585e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.36969e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20638e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.27238e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.31122e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.93659e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.56505e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.90742e-16
    FALSE [1] 1
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.9776e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.24558e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36849e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.59133e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.61859e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.83914e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16623e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07347e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.28978e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.87597e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.20636e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72338e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.19932e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.36201e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.23177e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.81813e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.87074e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.06441e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.67125e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.21286e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.29131e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.78263e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.97589e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.42584e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.31902e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.54104e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69166e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54482e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05634e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.66298e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.37709e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.15017e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.66196e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14385e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.12418e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52114e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.54626e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.56851e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.18087e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.50501e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.0282e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.0506e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.96582e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16586e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.86515e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.46848e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69316e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.66639e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.63489e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.75306e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13637e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63978e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.22124e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.66316e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67399e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.066e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.69345e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07586e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.84084e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.89999e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76254e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.93297e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.12431e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.21446e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.70949e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.42583e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.11768e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.06569e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07755e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.1963e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.45713e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63559e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.04654e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91353e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.61565e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72091e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.61939e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74605e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.44195e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.76029e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.0305e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.44113e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.14055e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.18945e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.0406e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.52235e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.81362e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.84309e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.0815e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.99299e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60493e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07195e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.0057e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.04467e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.99456e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25703e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50116e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.44583e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1124e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30467e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30977e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06568e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.50784e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03244e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.14115e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08003e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.9746e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.60974e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.86107e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4676e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07871e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11609e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.01999e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.73939e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15203e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.9454e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.70342e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.57796e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.5137e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47001e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16063e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.00417e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.91752e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.1201e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05202e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.2181e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.25193e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.80348e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16946e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2003e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98728e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19729e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.08326e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27765e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.4525e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.41601e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.875e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20781e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.53033e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.34546e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.9248e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.1069e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47718e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.77463e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.83583e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.21671e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07529e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.19644e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41446e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.51545e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.66175e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65519e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.06497e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.93864e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.28502e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27408e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.39852e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.94007e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.2149e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.91755e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.82954e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.29746e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.80733e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.43144e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30024e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01666e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.61409e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.67182e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.40675e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.89442e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.10722e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.88823e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.86217e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32025e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34136e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.93499e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.11665e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.04387e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98619e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.16813e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.96411e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.89205e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33891e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35932e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.94873e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.5667e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.26964e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73132e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.20885e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.19702e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.46094e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35648e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.37601e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00571e-17
    FALSE [1] 2
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08838e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20927e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.53286e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.62239e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.7266e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.46976e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.31084e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.77333e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.47342e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.14202e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.0612e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.64425e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.4253e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.49638e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02404e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.51466e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.18418e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.94866e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35915e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.36756e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.11577e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73519e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.94872e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1476e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07061e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09322e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.79192e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.80524e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34141e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.52067e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77526e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.51138e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.13753e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.93794e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.29827e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.66458e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16122e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.79686e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.39544e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16391e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.01457e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.66874e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.04115e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.02057e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69826e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.08239e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.69977e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28343e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59157e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.64598e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06192e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.93409e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.5493e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.80499e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.86177e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28081e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.08961e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33207e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.86134e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.04218e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35658e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.58041e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02106e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.91119e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.89044e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.67706e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.3537e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59742e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.50725e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.24976e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.15044e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.35053e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67502e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.29001e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34271e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19322e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.82095e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.56053e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.18196e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32496e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.10637e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.774e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4096e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.06904e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.51207e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.09541e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.05052e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32264e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.95187e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.08284e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38286e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28608e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.25787e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.94655e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.27333e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65257e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.36129e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33035e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49329e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2328e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.97695e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.70286e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.21223e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.56741e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16245e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.58549e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.4322e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08119e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4145e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.74753e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33231e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.29872e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.88126e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.66552e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15676e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.50333e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.37241e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.68312e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98666e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06734e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.82101e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.70401e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.94029e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.04231e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60435e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.35956e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.67512e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20758e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.09584e-17

    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.83584e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.04483e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77905e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.45607e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.27376e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01032e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.64145e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.2021e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15976e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.46642e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.137e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.9287e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.90814e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16662e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56964e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.27156e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00198e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.37369e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.44564e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.38998e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.3123e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.42127e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09932e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.62312e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.89155e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59941e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.44079e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54383e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.62692e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.27472e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.39898e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38114e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63746e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49029e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18329e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.37619e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4427e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.68021e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.10095e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73591e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28638e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.777e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50596e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.36761e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.16547e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38388e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.13549e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.93681e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.30756e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.34583e-17
    FALSE [1] 3
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.24885e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.51666e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.69236e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.03815e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26682e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.15351e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.37718e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33348e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27246e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.26622e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.62464e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.68535e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27814e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.08237e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.0846e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.30757e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16506e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58654e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.03327e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.50365e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.33844e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.564e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.29481e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19493e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.95097e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65505e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.15003e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.849e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.53701e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.95332e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.33279e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.31982e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54425e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14958e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.10999e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.02532e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.64913e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.37587e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74956e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.86085e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.31102e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01456e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.7126e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.70216e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.84991e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.91367e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65625e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.84683e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.55867e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.25363e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16619e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.64077e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.8979e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26076e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.93346e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.68454e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.06474e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18282e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01096e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.25108e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.19869e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.43592e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.69704e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20969e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.96203e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.96485e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20206e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.22261e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.01361e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.70998e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.48389e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.17637e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.13969e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05054e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17451e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.28853e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.43612e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07532e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02349e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1755e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25012e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.87206e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16754e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.86162e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.62966e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17036e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.58084e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08139e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35425e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.24426e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14546e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.29668e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.41322e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25337e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.00637e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03683e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.51676e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.8713e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.97083e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.84432e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.50447e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.10449e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.40822e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.00937e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.53481e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.45841e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.1216e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.78126e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.43699e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.29107e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36354e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.61644e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1126e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50985e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.04439e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.07494e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.47934e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.76156e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.26229e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.11582e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.44793e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.5675e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.72568e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06859e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.31301e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.56228e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00203e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.49944e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.00723e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73338e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.41343e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.92055e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.26342e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.31022e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.35831e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.04749e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.92651e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.184e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.11952e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06193e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.25088e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.52943e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.76338e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2582e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.32527e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.97222e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.10733e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.59578e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.73173e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.11354e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.03222e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.45858e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11301e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77353e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.01494e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.6114e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.64375e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38563e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.84148e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.53266e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58239e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.90778e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25347e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02919e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69745e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.4364e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.78442e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.21283e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.67958e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.57456e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.26641e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26328e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.31528e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.4919e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.44385e-17
    FALSE [1] 4
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.99254e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.43125e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48529e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74129e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.52445e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.64456e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74522e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.94996e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.68132e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.33847e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.91844e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.30885e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52558e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.04222e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.36411e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.3943e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.17687e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.60439e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.50056e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58958e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54062e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00055e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.37102e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2819e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.77657e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.9989e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.42135e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.84042e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.24861e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.42718e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.08607e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35709e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.0134e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.45912e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28303e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.31981e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.8975e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.58356e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.12678e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.31274e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.4863e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.04441e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02545e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.49599e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25523e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41986e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.90629e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.59806e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72253e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73344e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69877e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.97766e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60168e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.03565e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.97793e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.35738e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.79334e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03315e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.55875e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.87498e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.55299e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.7818e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.87581e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.29959e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.86052e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.09641e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.98954e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.0212e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.86834e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.62473e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.5564e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.1831e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.25364e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.35426e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.68054e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26124e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.15317e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.67963e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.72398e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.22822e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74061e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.94668e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.68013e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.23469e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.97763e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.93413e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.00776e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.00551e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.3973e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91795e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.31966e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.72221e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.3433e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.21543e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.04463e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.6979e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.40124e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33232e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.11846e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.29867e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.86263e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06196e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.81434e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4878e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41009e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.52584e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.88835e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.01244e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.06855e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54358e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.18254e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15702e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.32323e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52513e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.83273e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36561e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.74866e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16774e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05036e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13494e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.1282e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65379e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.11463e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.38394e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.90517e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.25539e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.23025e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41277e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76237e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.63402e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.7068e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.01475e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32805e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.21558e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.23002e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.14656e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09745e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56738e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.84848e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.24061e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.79285e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.0166e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.01378e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26535e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.8671e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.42454e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.72079e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.04811e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47054e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.50525e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50758e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01547e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74671e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77212e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.459e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.72164e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.85719e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.45612e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.93629e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00512e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.27476e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.5522e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76105e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.09946e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.37275e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18438e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36647e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.6503e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.00234e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.85624e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.3347e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.83841e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.80799e-17
    FALSE [1] 5
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.30776e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.08834e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50214e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.66286e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.70873e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.024e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.85941e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.33475e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.80823e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.22509e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.94885e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.99251e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.76709e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67808e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20133e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50865e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.43675e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.43525e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58803e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.51869e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.35068e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.1888e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.99343e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.79539e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.3317e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11492e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.49479e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.73511e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.99235e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.1446e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.05442e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.03134e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.33308e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77202e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91019e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.58026e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.58556e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.83693e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.54154e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.28233e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.87754e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.53408e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74609e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67187e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.68766e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.5804e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.82671e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.03533e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59718e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76189e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.79895e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91994e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.44107e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.45619e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.24224e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.21685e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.42694e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.9054e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.99523e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.27188e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00756e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.95333e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.75966e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76859e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.29823e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.48797e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36179e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.97722e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.74618e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.31262e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.47434e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.2108e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.39241e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.68352e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.82337e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.73867e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.0988e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.00733e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.62262e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.59716e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15749e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.54871e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17074e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.00214e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.45459e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35166e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.8739e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.87165e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.05004e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.70199e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.51452e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.19474e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.60038e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26681e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.21757e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.42191e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.75898e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.71948e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27498e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.5838e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11366e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.23323e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1614e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.33367e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.31985e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.90563e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.75822e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.78566e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.86077e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.85702e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.9756e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.78584e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63805e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.14751e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.42216e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.13763e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41591e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.76022e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.10243e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2851e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72526e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13635e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.42699e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.92009e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.36977e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.48416e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18108e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11259e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.88464e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.88616e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.90997e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20977e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.23426e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.22661e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.9306e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41744e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.78217e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.84285e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16123e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.66151e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.92979e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28809e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54007e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.43331e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.96307e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.00551e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.49162e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65427e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17063e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4682e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07467e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.69677e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.67689e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27546e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09915e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.57382e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.72106e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07403e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.38909e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.65789e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.72236e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.35429e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16883e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.82088e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.24447e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.80084e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.33818e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.05567e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.9632e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.06003e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.10577e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.81637e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15749e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.28322e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07827e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.43934e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69273e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.62164e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.75014e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.24153e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.17032e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03822e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.52723e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.11801e-16
    FALSE [1] 6
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.74521e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.72461e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.68103e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30124e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16633e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.16519e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.17008e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.10496e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.93213e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.99752e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.59835e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.68103e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20902e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.97625e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.26937e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.48715e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.46864e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.12019e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.0737e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02332e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60832e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.97939e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.05496e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91535e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.96214e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47296e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.57009e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.45251e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.57014e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.64532e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.74133e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.28982e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.28276e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.51827e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.9802e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.54326e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54986e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.58939e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.07396e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.51541e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.85312e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49535e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.62545e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.97644e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.79665e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.39487e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.47309e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.02668e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.66302e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.42064e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.32524e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.97281e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.26637e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60505e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.54889e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.39553e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.74634e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09774e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.3864e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.64635e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.00213e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.45949e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73669e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.59013e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01593e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.42351e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.55006e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.61247e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.52779e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.31202e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28217e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.16612e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.25289e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.74087e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.79919e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.34399e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48349e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.88549e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.42295e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.24499e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.05609e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.48567e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.46001e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.95302e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63966e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.84923e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35208e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.53793e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.18963e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.75282e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.15666e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41934e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.46223e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.22356e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.44205e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58272e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.24431e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03416e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.34901e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.22882e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.53626e-17

    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.80078e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.9526e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.15315e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.93158e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.18512e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.9618e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72556e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.53013e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16424e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.53207e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.23847e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.54685e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.64902e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.43628e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.05618e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.86607e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.79719e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.37874e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.84138e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48012e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.6827e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.93486e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.12245e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.38661e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.37383e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.30432e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00113e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.62799e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.54958e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.04505e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19507e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.11156e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48569e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.76419e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.86046e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.74433e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.40245e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.30745e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28849e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.58271e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50587e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.30791e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56277e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.66609e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.9624e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.37785e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36198e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.71185e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59402e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.73793e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30391e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.4207e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.71447e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.66257e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65264e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.46792e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.47948e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.32988e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.97867e-17

    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.89483e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.87846e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.28062e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00877e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.62219e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36634e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.57782e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.80881e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.79702e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.9527e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.70007e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.1749e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74065e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20221e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.85975e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.51058e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.83626e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.00768e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.93588e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.70916e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19507e-17
    FALSE [1] 7
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.03393e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.12077e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00376e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.70846e-20
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.13153e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38024e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.14202e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.21588e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.66567e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.08969e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.35264e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.63442e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.63802e-20
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.72088e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.33738e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.01152e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.21588e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01351e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16394e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.7897e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.50543e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54597e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.22541e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.48745e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.3315e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.43701e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.82042e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.63852e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.20894e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.84864e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74881e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56401e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.18498e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54566e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.0567e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.21774e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.94957e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.58724e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.87672e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65052e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52556e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.03463e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.88077e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.3731e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1852e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.5806e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.25418e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.32205e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.68815e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.39931e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.46782e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.3738e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.89408e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.54933e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.46199e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.12475e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.62818e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.10237e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.37885e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.57609e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.54862e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14894e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.18065e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20265e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.65196e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.14176e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05512e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20608e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59246e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.63523e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.98802e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.81152e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.59602e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.13473e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.40424e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.99095e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.0269e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28037e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.18988e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.21397e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.91407e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.51735e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.32268e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.83363e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.31388e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91001e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49231e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.69926e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30398e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07342e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.69097e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.1527e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47301e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.37434e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07857e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.43986e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.06595e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.86151e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.07472e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.40482e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76205e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.74193e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.67548e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.61042e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.06508e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.8082e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.16208e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.63269e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38262e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.78218e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03939e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.96333e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69444e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.49913e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.0617e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.84755e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.29828e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38851e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.0727e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48183e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.5897e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.10654e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.28912e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.68968e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49217e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.02347e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.65669e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54184e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13911e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05522e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.97206e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08765e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25688e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.56262e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.714e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.64783e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08049e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.39373e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69973e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.29441e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.81699e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17947e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.14033e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.42221e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16317e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.69415e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67726e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.27409e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.24214e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33431e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.7829e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.29631e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03037e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30207e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.41021e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.94617e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13732e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91176e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47156e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.1241e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.64645e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.48186e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33514e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.04095e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.88604e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.04569e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73191e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.58634e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28376e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52548e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.70619e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.82529e-20
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.96784e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08819e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69867e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.75258e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30049e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60053e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.1701e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32391e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.92e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.44196e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.2316e-18
    FALSE [1] 8
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.6337e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.75035e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.52256e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.21973e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.39058e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.40983e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00414e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.41977e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.79261e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.19514e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.66831e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.81797e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.61755e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.30452e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05142e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.9731e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41607e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.21771e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13224e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.37904e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.0228e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20824e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.33099e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.83894e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65668e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.51372e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.43958e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.7598e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.66753e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.19695e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.96179e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.99063e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.06107e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.5539e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.00085e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2185e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.87796e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.54788e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49884e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.83138e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.78516e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.62634e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.98241e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.067e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.19567e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38051e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.84221e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.34743e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.44447e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.54563e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20834e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.61599e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.0451e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.05217e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.48268e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48899e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01468e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.30488e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.48315e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.5851e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.00745e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32424e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.87531e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91206e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.7108e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56298e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.92587e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.30755e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.43681e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.53208e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.36717e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.99231e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65989e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02528e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.85633e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09276e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27134e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.88132e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.94243e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09719e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.88889e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.76092e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.93197e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.63399e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.11222e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.97781e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.24309e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.9033e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17403e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.68356e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.55346e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18936e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.87041e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.28373e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50577e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.64769e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.78442e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.23817e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.60135e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65478e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77457e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.52027e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.5637e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.61354e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.22514e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.59194e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17745e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.20506e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.49264e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.10518e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.16113e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08179e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.6171e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.14532e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34016e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.03814e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02529e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.8768e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.53988e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.17824e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.92474e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00678e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.82888e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.54811e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.21718e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.35012e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63389e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.22899e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.34536e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48798e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25453e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.63868e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.99492e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.32691e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65421e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.4847e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65211e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.40859e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.02549e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.39836e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.72355e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.34909e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.03351e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54709e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.83741e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.94645e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.31631e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32536e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.9107e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.11255e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52669e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00252e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27942e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.26359e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.54326e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.52358e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.67259e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16819e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.86871e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.6034e-20
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.07395e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.85266e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.28555e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73806e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.89492e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.31029e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.59134e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07844e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1448e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.93515e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19866e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.70089e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.91023e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15141e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.3895e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.15912e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.90517e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.08722e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.88099e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.4303e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.36944e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.81278e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.71392e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.46932e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.0446e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.61483e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.59793e-20
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.75238e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.07895e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.38904e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.43988e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.57234e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.29233e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.78658e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.2942e-20
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.70645e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.57373e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.64372e-17
    FALSE [1] 9
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.08498e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.84551e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.99715e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.10199e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.01994e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.91157e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00327e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.87401e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.91468e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.00673e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.31742e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69057e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.12839e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.50066e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.68894e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.3282e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.91285e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.82505e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.15684e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.40385e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28179e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19025e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01851e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.19764e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.1456e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2265e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.53642e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.47926e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.40664e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.83132e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.32784e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.9161e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.36657e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17231e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.95844e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.66496e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.19738e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.59046e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.57162e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.68743e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.50255e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.86833e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09225e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.00393e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02539e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.87239e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.72922e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.10377e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.93918e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03006e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.20403e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.24124e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77483e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.26141e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.95005e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.12661e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.67237e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.02905e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.38779e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.63268e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.97169e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.23723e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.18705e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.45661e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14127e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.27498e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.10923e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52537e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.54734e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41341e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.24303e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.35936e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.47675e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16102e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.0735e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.39515e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.82362e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.65544e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.75333e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52779e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.02283e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56152e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.82882e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.19381e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18873e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.91273e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98013e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.03309e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13472e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.14002e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.60148e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.66312e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07062e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.21172e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.51597e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.78847e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.92331e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.73846e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2049e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.19568e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.11455e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.88004e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.40866e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41077e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.00873e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11987e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32051e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.3052e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06099e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.143e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.80398e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.11493e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.44953e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.73034e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.63431e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.51707e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.79219e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.27585e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.80728e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.7391e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.64659e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.39045e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.49812e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.25921e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.83841e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17627e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.53757e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.00021e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.46465e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.56454e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.84165e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.30171e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72015e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.71703e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.54406e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.88332e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.59324e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48641e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.17753e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.69195e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.06391e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.67308e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.68748e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.53688e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08942e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.87731e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.98684e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.50032e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.5585e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.00988e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.19476e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.13256e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.23862e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28197e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.49058e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.99711e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.49385e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67676e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91021e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91146e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76373e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.2384e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98427e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.65372e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.23735e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.87286e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.33948e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.29499e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.9134e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.53094e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.56541e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.79763e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.73037e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.96721e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.7076e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.77242e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.30563e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.32352e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.72399e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.34975e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.01351e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.94003e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.13734e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.51035e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.42669e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.12542e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.2129e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.2178e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.68504e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.48086e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.61287e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60166e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.16236e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.09429e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.33429e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.8644e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.69229e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.60248e-16
    FALSE [1] 10
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.0834e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.30586e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41817e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.18875e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.65548e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.30193e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.75758e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.41504e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.52281e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.70077e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.03385e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.65296e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.33989e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.89689e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.98707e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.03575e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.86101e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.1449e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.94953e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.71083e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.87609e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.31889e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.02251e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.22985e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.31677e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.53483e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.04553e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.67477e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19795e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.17452e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.44961e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.49287e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.7281e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.72371e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.70092e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16139e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.16932e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.72586e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.13945e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.50913e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.13528e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.09791e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.28715e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.19099e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.66754e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.9704e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.08543e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.46228e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.15251e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.99557e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.52254e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91786e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.71227e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.88853e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.3159e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.31916e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.04397e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.20425e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.43762e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.92251e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.92612e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05105e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.11385e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.97043e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.72548e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.54482e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.24225e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.51873e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.21196e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.15369e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.32435e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.18678e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.0029e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20383e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.66599e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.20726e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.80518e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.83252e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.25719e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.06073e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.52603e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.01998e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.85685e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.0219e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.03823e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.16825e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07101e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.07389e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.35807e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.31199e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.02546e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.18475e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.13189e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.11091e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.4357e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.78808e-19
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.89322e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.9883e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.38317e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.47518e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.13788e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.64455e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.35875e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.90523e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.31788e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.05544e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.98896e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.99854e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.14457e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.93108e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.83849e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.80692e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.71671e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.7556e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.17939e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.88321e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07377e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.6148e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.37994e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.28614e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.70601e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.47146e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.07143e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.45259e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.19013e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.08182e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.42343e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.02209e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.84113e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.945e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.20338e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.11607e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.50405e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.91247e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.26728e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.22318e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.41742e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.29691e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.91497e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.5713e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.88109e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.9735e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.14249e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.48419e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.80654e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.86954e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.55658e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.16829e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.02984e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.23134e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.99683e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.14754e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.04197e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.72249e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.13054e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.41855e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 9.73442e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.52807e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.10352e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.77176e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.90252e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.09809e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.75995e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.27626e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 5.40304e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.2898e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.20511e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.66464e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.01738e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.36326e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.56888e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.45591e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.58905e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.71989e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.55042e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.95739e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.30979e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.01394e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.13834e-18
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.60174e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.42528e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.07922e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.83653e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.91945e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.78813e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.13952e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.20497e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.21622e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.24488e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.6901e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.28788e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.00216e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.65753e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 7.80787e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 8.44559e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.4903e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.01374e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 2.49142e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.76973e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.15069e-16
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 6.21664e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 4.78317e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 3.05749e-17
    FALSE Error in solve.default(AUX1 + 2 * Sigmalambda) : 
    FALSE   sistema es computacionalmente singular: número de condición recíproco = 1.90742e-16

    FALSE    user  system elapsed 
    FALSE 703.393   0.087 703.443

The results obtained are the following:

``` r
fit.rob$nknots
```

    ## [1] 0

``` r
fit.rob$coef.const
```

    ## (Intercept) 
    ##  -0.1552445

``` r
fit.rob$coef.lin
```

    ##  [1]  3.266113e-01 -6.285045e-04 -4.492650e-01 -1.443667e-01  3.136281e-02
    ##  [6]  9.270391e-04 -1.576092e-04 -1.527110e-08 -6.602658e-14  1.104379e-01

``` r
fit.rob$la1
```

    ## [1] 0.02

``` r
fit.rob$la2
```

    ## [1] 1.3

Are there covariables irrelevant for the model?

``` r
fit.rob$is.zero
```

    ##  [1] FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
    ## [13]  TRUE

It can be appreciated that only five covariates are considered important
for the model for the robust approach.

Let’s see if the proposal identifies any large residuals:
![](README_files/figure-markdown_github/residuals-1.png)

    ## [1] 17

17 observations were identified by the boxplot. These observations
corresponds to observations

``` r
in.ro <- (1:length(res))[ res %in% aa$out ]
in.ro
```

    ##  [1]  28  35  39  40 137 148 163 168 178 182 208 219 223 262 263 270 299

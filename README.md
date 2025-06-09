In this package, a robust estimation and variable selection procedure for partially linear additive models is performed using the real dataset from nutritional epidemiology.

Let's first install package <code>rplam-vs</code>

``` r
library(devtools)
install_github("alemermartinez/rplam-vs")
library(rplam-vs)
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


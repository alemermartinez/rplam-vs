######################################
#--- Plasma beta-carotene example ---#
######################################

#Read the data
datos <- read.table("Plasma_Retinol.txt", sep="\t")

#Overview of variables
str(datos)

#Standarizing robustly the continuos variables

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

#Response variable
betaplasma <- (datos$V13-median(datos$V13))/mad(datos$V13)

#Comparison
pairs(datos)

#Model
y <- betaplasma
Z <- cbind(sex,smok1,smok2,quetelet,vit1,vit2,calories,
           fat,alcohol,betadiet)
X <- cbind(age,cholesterol,fiber)

#Extreme outlier in alcohol:
boxplot(alcohol)

#Read the script with functions
source("R/rplam-vs-fn.R")

#Computing the robust estimator with the complete dataset

#- Grid used for the 10 covariates of the linear part
grilla.chica <- c(0, 0.05, 0.1)
grid.lambda <- expand.grid(grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica)
n.grid <- length(grilla.chica)

#- Fit the robust estimator 
#- IMPORTANT: Since the large number of covariates, this is a long procedure.
#- 6660 segs aprox in a Intel® Core™ i7-10700 CPU @ 2.90GHz × 16
set.seed(123)
{system.time(
  fit.full <- plam.rob.vs.betas.lambdas(y= y, Z= Z, X= X, grid.lambda=grid.lambda)
)}

#- Adjusted values
y.adjusted <- fit.full$prediction

#- Residuals
res <- y - y.adjusted

#- Plots of the residuals
plot(res)
hist(res)

#- Boxplot
aa <- boxplot(res)

#- Outliers detection
in.ro <- (1:length(res))[ res %in% aa$out ]
in.ro

#- Plot of the adjusted values vs residuals
plot(y.adjusted, res, xlab="Adjusted values", ylab="Residuals")




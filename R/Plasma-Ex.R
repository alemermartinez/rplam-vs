#Ejemplo de Plasma Retinol


datos <- read.table("Plasma_Retinol.txt", sep="\t")
str(datos)

#Si bien lo probé sin la estandarización (que sí selecciona más variables pero los errores, obviamente, son más grandes)
#lo "correcto" sería con la estandarización sólo para poder compararlo con lo hecho en los otros papers.

age <- scale(datos$V1) #datos$V1 #scale(datos$V1)
sex <- datos$V2-1
smokstat <- datos$V3
smok1 <- as.numeric(datos$V3==2)
smok2 <- as.numeric(datos$V3==3)
quetelet <- scale(datos$V4) #datos$V4 #scale(datos$V4)
vituse <- datos$V5
vit1 <- as.numeric(datos$V5==1)
vit2 <- as.numeric(datos$V5==2)
calories <- scale(datos$V6) #datos$V6 #scale(datos$V6)
fat <- scale(datos$V7) #datos$V7 #scale(datos$V7)
fiber <- scale(datos$V8) #datos$V8 #scale(datos$V8)
alcohol <- scale(datos$V9) #datos$V9 #scale(datos$V9)
cholesterol <- scale(datos$V10) #datos$V10 #scale(datos$V10)
betadiet <- scale(datos$V11) #datos$V11 #scale(datos$V11)
retdiet <- scale(datos$V12) #datos$V12 #scale(datos$V12)
retplasma <- scale(datos$V14) #datos$V14 #scale(datos$V14)

#Response variable
betaplasma <- scale(datos$V13) #datos$V13 #scale(datos$V13)

pairs(datos)


#Mismo modelo que Guo2013
y <- betaplasma
Z <- cbind(sex,smok1,smok2,quetelet,vit1,vit2,calories,
           fat,alcohol,betadiet)
X <- cbind(age,cholesterol,fiber)

#-- Análisis de outliers --#
#Alcohol:
boxplot(datos$V9)
#Hay un outlier extremo

#Edad
boxplot(datos$V1) #No hay outliers

#Quetelet
boxplot(datos$V4) #Hay outliers pero todos pegados

#Calories
boxplot(datos$V6) #También tiene un outlier bastante lejano al resto de los outliers

#Fat
boxplot(datos$V7) #Hay outliers pero todos bastante cercanos entre sí

#Fiber
boxplot(datos$V8) #Hay outliers pero todos bastante cercanos entre sí

#Cholesterol
boxplot(datos$V10) #Hay outliers pero todos bastante cercanos entre sí

#Betadiet
boxplot(datos$V11) #Hay outliers pero todos bastante cercanos entre sí

### Estimador

degree.spline <- 3

source("R/rplam-vs-fn.R")
library(fda)

set.seed(123)

#Grilla
grilla.chica <- c(0,0.05,0.1) #c(0, 0.001, 0.01, 0.05, 0.1, 0.2)
grid.lambda <- expand.grid(grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica)
dim(grid.lambda)

system.time(
  fit.rob <- plam.rob.vs.betas.lambdas(y, Z, X, degree.spline=degree.spline, grid.lambda)
)#Con dos elementos en la grilla tarda 133.179

fit.rob$nknots
fit.rob$coef.const
mean(y)
fit.rob$coef.lin
fit.rob$lambdas
fit.rob$is.zero
median( abs(y-fit.rob$prediction))

aa <- boxplot(y-fit.rob$prediction)
length(aa$out)

#Gráficos con residuos parciales
for(i in 1:10){
  p.res <- y-fit.rob$coef.const-Z[,-i]%*%as.matrix(fit.rob$coef.lin[-i])-rowSums(fit.rob$g.matrix)
  plot(Z[,i],p.res)
  ord <- order(Z[,i])
  lines(Z[ord,i],Z[ord,i]%*%as.matrix(fit.rob$coef.lin[i]),col='blue',lwd=2)
}

for(i in 1:3){
  p.res <- y-fit.rob$coef.const-Z%*%as.matrix(fit.rob$coef.lin)-rowSums(fit.rob$g.matrix[,-i])
  plot(X[,i],p.res)
  ord <- order(X[,i])
  lines(X[ord,i],fit.rob$g.matrix[ord,i],col='blue',lwd=2)
}


for(i in 1:3){
  ord <- order(X[,i])
  plot(X[ord,i],fit.rob$g.matrix[ord,i],type="l",col='blue',lwd=2)
}


############################################
# Acá comienza la iteración con el robusto #
############################################

source("R/rplam-vs-fn.R")
grilla.chica <- c(0, 0.1)
grid.lambda <- expand.grid(grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica,grilla.chica)

n <- length(y)
M <- 10 #Debería ser 100 (o 20)
contador.max <- 1000
mar <- rep(NA,M)
mape <- rep(NA,M)
sizes <- rep(NA,M)
nknots.muestras <- rep(NA,M)
lambdas.muestras <- matrix(NA,M,10)
cant.out <- rep(NA,M)

set.seed(17)
submuestras <- matrix(NA,M,100)
is.zero <- matrix(NA,M,13)
g1 <- matrix(NA,M,215)
g2 <- matrix(NA,M,215)
g3 <- matrix(NA,M,215)

#Lo que sigue es lento!
runs <- 0
contador <- 0
system.time({
while( (runs<M) & (contador<contador.max)) {
  contador <- contador+1
  subsample <- sample(1:n,100,replace=F)
  if( (range(X[-subsample,1])[1]<=range(X[subsample,1])[1]) &
      (range(X[-subsample,1])[2]>=range(X[subsample,1])[2]) &
      (range(X[-subsample,2])[1]<=range(X[subsample,2])[1]) &
      (range(X[-subsample,2])[2]>=range(X[subsample,2])[2]) &
      (range(X[-subsample,3])[1]<=range(X[subsample,3])[1]) &
      (range(X[-subsample,3])[2]>=range(X[subsample,3])[2]) ){


      print(contador)
      runs <- runs+1
      print(runs)

      y.new <- y[-subsample]
      Z.new <- Z[-subsample,]
      X.new <- X[-subsample,]
      np.point.new <- X[subsample,]
      fit.full <- plam.rob.vs.betas.lambdas(y= y.new, Z= Z.new, X= X.new,
                                             np.point = np.point.new, grid.lambda=grid.lambda)
      #if (class(fit.full)[1] != "try-error") {



        is.zero[runs,] <- fit.full$is.zero
        nknots.muestras[runs] <- fit.full$nknots
        lambdas.muestras[runs,] <- fit.full$lambdas
        print(lambdas.muestras[runs,])

        submuestras[runs,] <- subsample
        g1[runs,] <- fit.full$g.matrix[,1]
        g2[runs,] <- fit.full$g.matrix[,2]
        g3[runs,] <- fit.full$g.matrix[,3]

        predicciones <- fit.full$coef.const + as.vector(fit.full$coef.lin%*%t(Z[subsample,])) + rowSums(fit.full$np.prediction)
        errores <- y[subsample]-predicciones
        mape[runs] <- median(abs(errores))

        aa <- boxplot(errores)
        cant.out[runs] <- length(aa$out)

        residuos <- y[-subsample]-fit.full$prediction
        mar[runs] <- median(abs(residuos))
        sizes[runs] <- sum(fit.full$is.zero==FALSE)
      #}
  }
}
}) #Tarda 287.73 segs=4.79 mins


#Cuáles son cero:
is.zero.rob <- colSums(is.zero)

#Average size:
sizes
av.size.rob <- mean(sizes)

#MAR:
mar
mar.rob <- mean(mar)

#MAPE:
mape
mape.rob <- mean(mape)

#Selected frequency for the 10 variables as 0
col.is.zero.rob <- colMeans(is.zero)

#Outliers
mean(cant.out)
median(cant.out)

#Grafico outliers
nombre= "Cantidad-Outliers-100Muestras.pdf"
pdf(nombre, bg="transparent")
hist(cant.out, xlab="cantidad de outliers", main="")
dev.off()


################################
# Ver qué queda con el clásico #
################################

grid.la1 <- c(0,0.001,0.01,0.05, 0.1, 0.2) #c(0,0.001,0.01,0.05, 0.1, 0.2,1,5,20) #c(0,0.1,1,10,100,1000) #c(0,0.001,0.01,0.05, 0.1, 0.2) #Con estandarizacion
grid.la2 <- c(0.01,0.05, 0.1, 0.2, 0.5, 1, 2) #c(0,0.1,1,10,100,1000) #10 #c(0,0.001,0.01,0.05, 0.1, 0.2) #Con estandarizacion

n <- length(y)
M <- 100 #Debería ser 20
contador.max <- 1000
mar <- rep(NA,M)
mape <- rep(NA,M)
sizes <- rep(NA,M)
nknots.muestras <- rep(NA,M)
la1.muestras <- rep(NA,M)
la2.muestras <- rep(NA,M)

set.seed(17)
submuestras <- matrix(NA,M,100)
is.zero <- matrix(NA,M,13)
g1 <- matrix(NA,M,215)
g2 <- matrix(NA,M,215)
g3 <- matrix(NA,M,215)

runs <- 0
contador <- 0
system.time({
  while( (runs<M) & (contador<contador.max)) {
    contador <- contador+1
    subsample <- sample(1:n,100,replace=F)
    if( (range(X[-subsample,1])[1]<=range(X[subsample,1])[1]) &
        (range(X[-subsample,1])[2]>=range(X[subsample,1])[2]) &
        (range(X[-subsample,2])[1]<=range(X[subsample,2])[1]) &
        (range(X[-subsample,2])[2]>=range(X[subsample,2])[2]) &
        (range(X[-subsample,3])[1]<=range(X[subsample,3])[1]) &
        (range(X[-subsample,3])[2]>=range(X[subsample,3])[2]) ){


      print(contador)
      runs <- runs+1
      print(runs)

      fit.full <- plam.cl.vs(y= y[-subsample], Z= Z[-subsample,], X= X[-subsample,],
                              np.point = X[subsample,], grid.la1=grid.la1,
                              grid.la2=grid.la2, k.malos.max = 2)
      #if (class(fit.full)[1] != "try-error") {



      is.zero[runs,] <- fit.full$is.zero
      nknots.muestras[runs] <- fit.full$nknots
      la1.muestras[runs] <- fit.full$la1
      la2.muestras[runs] <- fit.full$la2

      submuestras[runs,] <- subsample
      g1[runs,] <- fit.full$g.matrix[,1]
      g2[runs,] <- fit.full$g.matrix[,2]
      g3[runs,] <- fit.full$g.matrix[,3]

      predicciones <- fit.full$coef.const + as.vector(fit.full$coef.lin%*%t(Z[subsample,])) + rowSums(fit.full$np.prediction)
      errores <- y[subsample]-predicciones
      mape[runs] <- median(abs(errores))

      residuos <- y[-subsample]-fit.full$prediction
      mar[runs] <- median(abs(residuos))
      sizes[runs] <- sum(fit.full$is.zero==FALSE)
      #}
    }
  }
}) #Tarda 38.987 segs

ncl <- nknots.muestras
la1cl <- la1.muestras
la2cl <- la2.muestras

#Cuáles son cero:
is.zero.cl <- colSums(is.zero)

#Average size:
sizes
av.size.cl <- mean(sizes)

#MAR:
mar
mar.cl <- mean(mar)

#MAPE:
mape
mape.cl <- mean(mape)

#Selected frequency for the 13 variables as 0
col.is.zero.cl <- colMeans(is.zero)




#############################
# Ahora con el full robusto #
#############################

n <- length(y)
M <- 100 #Debería ser 20
contador.max <- 1000
mar <- rep(NA,M)
mape <- rep(NA,M)
sizes <- rep(NA,M)
nknots.muestras <- rep(NA,M)

set.seed(17)
submuestras <- matrix(NA,M,100)
g1 <- matrix(NA,M,215)
g2 <- matrix(NA,M,215)
g3 <- matrix(NA,M,215)

runs <- 0
contador <- 0
system.time({
  while( (runs<M) & (contador<contador.max)) {
    contador <- contador+1
    subsample <- sample(1:n,100,replace=F)
    if( (range(X[-subsample,1])[1]<=range(X[subsample,1])[1]) &
        (range(X[-subsample,1])[2]>=range(X[subsample,1])[2]) &
        (range(X[-subsample,2])[1]<=range(X[subsample,2])[1]) &
        (range(X[-subsample,2])[2]>=range(X[subsample,2])[2]) &
        (range(X[-subsample,3])[1]<=range(X[subsample,3])[1]) &
        (range(X[-subsample,3])[2]>=range(X[subsample,3])[2]) ){


      print(contador)
      runs <- runs+1
      print(runs)

      fit.full <- plam.rob(y= y[-subsample], Z= Z[-subsample,], X= X[-subsample,],
                              np.point = X[subsample,])
      #if (class(fit.full)[1] != "try-error") {

      nknots.muestras[runs] <- fit.full$nknots

      submuestras[runs,] <- subsample
      g1[runs,] <- fit.full$g.matrix[,1]
      g2[runs,] <- fit.full$g.matrix[,2]
      g3[runs,] <- fit.full$g.matrix[,3]

      predicciones <- fit.full$coef.const + as.vector(fit.full$coef.lin%*%t(Z[subsample,])) + rowSums(fit.full$np.prediction)
      errores <- y[subsample]-predicciones
      mape[runs] <- median(abs(errores))

      residuos <- y[-subsample]-fit.full$prediction
      mar[runs] <- median(abs(residuos))
      #}
    }
  }
}) #Tarda 28.955

nrob.full <- nknots.muestras

#MAR:
mar
mar.rob.full <- mean(mar)

#MAPE:
mape
mape.rob.full <- mean(mape)



#############################
# Ahora con el full clasico #
#############################

n <- length(y)
M <- 100 #Debería ser 20
contador.max <- 1000
mar <- rep(NA,M)
mape <- rep(NA,M)
sizes <- rep(NA,M)
nknots.muestras <- rep(NA,M)

set.seed(17)
submuestras <- matrix(NA,M,100)
g1 <- matrix(NA,M,215)
g2 <- matrix(NA,M,215)
g3 <- matrix(NA,M,215)

runs <- 0
contador <- 0
system.time({
  while( (runs<M) & (contador<contador.max)) {
    contador <- contador+1
    subsample <- sample(1:n,100,replace=F)
    if( (range(X[-subsample,1])[1]<=range(X[subsample,1])[1]) &
        (range(X[-subsample,1])[2]>=range(X[subsample,1])[2]) &
        (range(X[-subsample,2])[1]<=range(X[subsample,2])[1]) &
        (range(X[-subsample,2])[2]>=range(X[subsample,2])[2]) &
        (range(X[-subsample,3])[1]<=range(X[subsample,3])[1]) &
        (range(X[-subsample,3])[2]>=range(X[subsample,3])[2]) ){


      print(contador)
      runs <- runs+1
      print(runs)

      fit.full <- plam.cl(y= y[-subsample], Z= Z[-subsample,], X= X[-subsample,],
                           np.point = X[subsample,])
      #if (class(fit.full)[1] != "try-error") {

      nknots.muestras[runs] <- fit.full$nknots

      submuestras[runs,] <- subsample
      g1[runs,] <- fit.full$g.matrix[,1]
      g2[runs,] <- fit.full$g.matrix[,2]
      g3[runs,] <- fit.full$g.matrix[,3]

      predicciones <- fit.full$coef.const + as.vector(fit.full$coef.lin%*%t(Z[subsample,])) + rowSums(fit.full$np.prediction)
      errores <- y[subsample]-predicciones
      mape[runs] <- median(abs(errores))

      residuos <- y[-subsample]-fit.full$prediction
      mar[runs] <- median(abs(residuos))
      #}
    }
  }
}) #Tarda 0.718

ncl.full <- nknots.muestras

#MAR:
mar
mar.cl.full <- mean(mar)

#MAPE:
mape
mape.cl.full <- mean(mape)

AA <- cbind(rbind(mar.rob, mar.cl, mar.rob.full, mar.cl.full), rbind(mape.rob, mape.cl, mape.rob.full, mape.cl.full),  c(av.size.rob, av.size.cl, 13,13))
colnames(AA) <- c("MAR", "MAPE", "AV.SIZE")
AA

#Para la semilla set.seed(24) y M=20
#> AA
#MAR      MAPE AV.SIZE
#mar.rob      0.3246898 0.3351997    4.95
#mar.cl       0.4157613 0.4247257    6.00
#mar.rob.full 0.2608975 0.2977610   13.00
#mar.cl.full  0.3698180 0.4062411   13.00

#Para la semilla set.seed(17) y M=20
#> AA
#MAR      MAPE AV.SIZE
#mar.rob      0.3127855 0.3254017    4.80
#mar.cl       0.4108202 0.4199546    6.45
#mar.rob.full 0.2630941 0.2996844   13.00
#mar.cl.full  0.3672353 0.4072200   13.00

#Para la semilla set.seed(17) y M=100
#> AA
#MAR      MAPE AV.SIZE
#mar.rob      0.3181265 0.3297196    4.33
#mar.cl       0.4108586 0.4268663    6.18
#mar.rob.full 0.2620123 0.3088842   13.00
#mar.cl.full  0.3669361 0.4048766   13.00


#> round(AA,3)
#MAR  MAPE AV.SIZE
#mar.rob      0.318 0.330    4.33
#mar.cl       0.411 0.427    6.18
#mar.rob.full 0.262 0.309   13.00
#mar.cl.full  0.367 0.405   13.00

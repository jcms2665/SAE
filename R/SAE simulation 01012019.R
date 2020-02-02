# Simulation tools for SAE
# Kreutzmann et al
# 03012020


options(warn=-1)
revisar<- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
  
}

p<-c("foreign","ggplot2","sae","survey","survey","questionr")
revisar(p)
rm(list=ls()) 





# III. Métodos directos para la desagregación de datos de pobreza

# Individuos que residen en las provincias
# prov
# provlab   ->   Nombre de la provincia

#Datos del ingreso
data(incomedata)
attach(incomedata)
data(sizeprov)

table(incomedata$provlab)
table(incomedata$prov)

n<-dim(incomedata)[1]                   # Tamaño muestral total
D<-length(unique(incomedata$prov))      # Numero de provincias (Ã¡reas o dominios)
nd<-as.vector(table(incomedata$prov))   # n Tamanos muestrales de las provincias
Nd<-sizeprov$Nd                         # N Tamanos poblacionales de las provincias

#El estimador es el umbral de pobreza, que
z<-6557.143
poor<-numeric(n)
poor[incomedata$income<z]<-1

#prueba

povinc.dir.res<-direct(y=poor,dom=prov,sweight=weight,domsize=sizeprov[,-1])
print(povinc.dir.res,row.names=F)

povinc.dir<-povinc.dir.res$Direct
povinc.dir.cv<-povinc.dir.res$CV
sum(povinc.dir.cv>20)


# Estimadores GREG

rm(list=ls())

library(sae) 
data(incomedata) 
attach(incomedata) 
data(sizeprov)

data(sizeprovage) 
data(sizeprovedu) 
data(sizeprovlab)

Nd<-sizeprov[,3]
Ndage<-as.matrix(sizeprovage[,-c(1,2)]) 
Ndedu<-as.matrix(sizeprovedu[,-c(1,2)]) 
Ndlab<-as.matrix(sizeprovlab[,-c(1,2)])
Pdage<-Ndage/Nd 
Pdedu<-Ndedu/Nd 
Pdlab<-Ndlab/Nd

a1<-Pdage[,3:5]
a2<-Pdedu[,c(2,4)]
a3<-data.frame(Pdlab[,2])
  


X<-cbind(const=rep(1,D),Pdage[,3:5],Pdedu[,c(2,4)],Pdlab[,2]) 
Xtot<-model.matrix(poor~age3+age4+age5+educ1+educ3+labor1)

provl<-unique(prov) # Índice de cada provincia 
p<-dim(Xtot)[2]
# Número de variables auxiliares
betad<-matrix(0,nr=D,nc=p) # Matriz con los coeficientes de regresión # para cada provincia (en filas)
Xd.est<-matrix(0,nr=D,nc=p) # Matriz de estimadores directos de las medias # de las variables auxiliares para cada provincia
povinc.greg<-numeric(D) # Vector de estimadores GREG en las provincias
povinc.greg.var<-numeric(D) # Vector con las varianzas estimadas # bajo el diseño de los estimadores GREG
for (d in 1:D){
  Xd<-Xtot[prov==provl[d],] # Valores de las variables auxiliares # para los indiv. de la provincia
  wd<-weight[prov==provl[d]] # Pesos muestrales para los individuos # de la provincia
  yd<-poor[prov==provl[d]] # Valores de la variable de interés # para los individuos de la provincia
  # Ajustamos la regresión para la provincia, con los pesos muestrales 
  betad[d,]<-coef(summary(lm(yd~-1+Xd, weights=wd)))[,1] # Estimadores directos de las medias de las variables auxiliares en la provincia
  Xd.est[d,]<-colSums(diag(wd)%*%Xd)/Nd[d]
  # Estimador GREG de la incidencia de pobreza para la provincia 
  povinc.greg[d]<-povinc.dir [d]+sum((X[d,]-Xd.est[d,])*betad[d,])
  # Varianza estimada bajo el diseño del estimador # GREG de la incidencia de pobreza 
  gd<-matrix(1/Nd[d]+ +(X[d,]-Xd.est[d,])%*%solve(t(Xd)%*%diag(wd)%*%Xd)%*%t(Xd),nr=nd[d]) ed<-yd-Xd%*%as.matrix(betad[d,],nr=p) povinc.greg.var[d]<-sum(wd*(wd-1)*(gd*ed)^2)
}
#

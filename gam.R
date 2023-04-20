# packages
install.packages("data.table")    # excel
install.packages("corrplot")      # Correlations
install.packages("dplyr")
install.packages("gam")
install.packages("lmtest")

# clearing environment 
rm(list = ls())
graphics.off()
# ctrl + L => to clear console

# Add Libraries in order to use them 
library(data.table)
library("corrplot")
library("ggplot2")                     
library("GGally")
library(pastecs)
library(psych)
library(dplyr)
library(gam)
library(lmtest)



data <- fread("uni.csv")
data <- data.frame(data[,2:9])
colnames(data)[1] = "GRE"
colnames(data)[2] = "TOEFL"
colnames(data)[3] = "UniRatings"
colnames(data)[8] = "Admit"
summary(data)

set.seed(1)

detect_outlier <- function(x) {
  Quantile1 <- quantile(x, probs=.25)
  Quantile3 <- quantile(x, probs=.75)
  IQR = Quantile3 - Quantile1
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}


remove_outlier <- function(dataframe) {
  columns <- names(dataframe)
  for (col in columns) {
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  print("Remove outliers")
  print(dataframe)
}

for(i in 1:2){
  data <- remove_outlier(data)
}

data <- data %>% 
  mutate_at(vars(Research, LOR, SOP, UniRatings), as.factor)

# RMSE di Train e di test

RMSE_function <- function(m) {
  yhat <- predict(m,data = data[train])
  RMSE_train <- sqrt(mean((data$Admit[train]-yhat)^2))
  
  
  yhat <- predict(m, newdata = data[-train,])
  RMSE_test <- sqrt(mean((data$Admit[-train]-yhat)^2))
  
  m$RMSE_train <- RMSE_train
  m$RMSE_test <- RMSE_test
  return(m)
}


train = sample(dim(data)[1], dim(data)[1]*0.7, replace = FALSE)

#modello lineare

model.lm <- lm(Admit~.-SOP-UniRatings,data=data, subset=train)
summary(model.lm)
RMSE_function(model.lm)

# Creare un modello GAM completo
model.gam0 <- gam(Admit~s(TOEFL)+s(CGPA)+s(GRE)+Research+SOP+LOR+UniRatings,data=data,subset = train)
model.gam0 <- RMSE_function(model.gam0)
summary(model.gam0)
print(model.gam0$RMSE_train)
print(model.gam0$RMSE_test)
par(mfrow=c(3,3))
plot(model.gam0,se=TRUE,rug=TRUE)
bptest(model.gam0)

#Rimossi regressori non significativi
model.gam1 <- gam(Admit ~ s(TOEFL,5) + s(CGPA,5) + s(GRE,5) +Research + LOR,data=data,subset = train)
par(mfrow=c(2,3))
plot(model.gam1,se=TRUE,rug=TRUE)
summary(model.gam1)
RMSE_function(model.gam1)
bptest(model.gam1)



# Ricerca del numero di nodi migliore
models = list()
for(i in 1:5){
  for(j in 1:5){
    for(l in 1:5){
      model <- gam(Admit~s(TOEFL,i)+s(CGPA,j)+s(GRE,l)+Research,data=data,subset = train)
      model <- RMSE_function(model)
      model$index <- c(i,j,l)
      models <- append(models,list(model))
    }
  }
}

# ritorna il valore minimo del RMSE di test e il valore dei dt per ogni regressori continuo
minRMSE <- function(){
  min <- models[[1]]$RMSE_test
  index <- 0
  for(i in 1:length(models)){
    if(models[[i]]$RMSE_test < min){
      min <- models[[i]]$RMSE_test
      index <- i
    }
  }
  print(min)
  print(models[[index]]$index)
}

minRMSE()


#Final model
model.finalgam <- gam(Admit~s(TOEFL,1)+s(CGPA,5)+s(GRE,1)+Research+LOR,data=data,subset = train)
par(mfrow=c(2,3))
plot(model.finalgam,se=TRUE,rug=TRUE)
summary(model.finalgam)
RMSE_function(model.finalgam)






# Analisi dei residui

shapiro.test(model.finalgam$residuals)
bptest(model.finalgam)

par(mfrow = c(3,3))
hist(model.finalgam$residuals,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(model.finalgam$residuals, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(model.finalgam$residuals),digits = 4),
                                      "- var:", round(var(model.finalgam$residuals),digits = 4)))
abline(c(0,0),c(0,length(model.finalgam$residuals)), col= "red", lwd = 2)

boxplot(model.finalgam$residuals)$out

qqnorm(model.finalgam$residuals, main='Residuals')
qqline(model.finalgam$residuals)

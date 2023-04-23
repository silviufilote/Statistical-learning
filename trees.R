# packages
install.packages("MASS")
install.packages("ISLR2")
install.packages("data.table")    # excel
install.packages("corrplot")      # Correlations
install.packages("ggplot2")       # Pairs Plots Correlations
install.packages("GGally")
install.packages("pastecs")
install.packages("psych")
install.packages("dplyr")
install.packages("performance")
install.packages("lmtest")

# clearing environment 
rm(list = ls())
graphics.off()

# Add Libraries in order to use them 
library(MASS)
library(ISLR2)
library(data.table)
library("corrplot")
library("ggplot2")                     
library("GGally")
library(pastecs)
library(psych)
library(dplyr)
library(performance)
library(lmtest)


# importing the dataset
data <- read.csv("C:\\Users\\fsilv\\Dropbox\\UNI\\stistics\\1.csv")
data <- data.frame(data[,2:9])
colnames(data)[1] = "GRE"
colnames(data)[2] = "TOEFL"
colnames(data)[3] = "UniRatings"
colnames(data)[8] = "Admit"


# removing outliers
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

# generating the dummy variables
data <- data %>% 
  mutate_at(vars(Research, LOR, SOP, UniRatings), as.factor)


# Train dataset
set.seed(2)
train = sample(dim(data)[1], dim(data)[1]*0.7, replace = FALSE)


################################################################################
#########################################  Tree
################################################################################

install.packages("tree")
library(tree)

tree_model <- tree( Admit ~ . , data, subset = train, split = "gini")
summary(tree_model)
plot(tree_model)
text(tree_model,pretty = 0)


# cross-validation for pruning
model_tree_cv <- cv.tree(tree_model,  FUN = prune.tree)
best = min(model_tree_cv$size[model_tree_cv$dev == min(model_tree_cv$dev)])
prune <- prune.tree(tree_model, best = best)
summary(prune)


# plotting tree
par(mfrow = c(1,2))
# dev = cross-validation error rate
# size = number of terminal nodes of each tree -> complexity
plot(model_tree_cv$size, model_tree_cv$dev, ylab = "cross-validation total rate", xlab = "Complexity = foglie", main ="Cross-validation total leaves")
abline(v = 8, col = "blue", lty = "dashed")
text(8, 15, "Benchmark", adj = c(0.5, -0.5))
plot(prune, main="Main title",
     xlab="X axis title",
     ylab="Y axis title",
     sub="Sub-title")
text(prune, pretty = 0, main = "Tree pruned")


# validation and MSE
yhat <- predict(tree_model, newdata = data[-train,])
res = data$Admit[-train] - yhat
test_RMSE_tree = round(sqrt(mean((res)^2)), digits = 3)



# residual analysis
bptest(tree_model_train) # non sono omoschedastici
shapiro.test(res) # non sono normali

par(mfrow = c(1,4))
hist(res,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(res, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(res),digits = 4),
                                      "- var:", round(var(res),digits = 4)))
abline(c(0,0),c(0,length(res)), col= "red", lwd = 2)

boxplot(res, ylab = "Residuals", main = "Outliers")$out

qqnorm(res, main='Residuals')
qqline(res)


################################################################################
#########################################  Bagging model
################################################################################

install.packages("randomForest")
library(randomForest)

set.seed(1)

bagg_test <- randomForest(Admit ~ . ,data = data, subset = train,
                           mtry = ncol(data)-1, importance = TRUE, replace = TRUE)
bagg_test
plot(bagg_test, main = "Bagged Trees: Error vs Number of Trees")
abline(v = 100, col = "blue", lty = "dashed")



# fare crossvalidazione sul numero di alberi
rmse_cv_bag = rep(0,10)
for (i in 1:10){
  bagg_model_cv <- randomForest(Admit ~ . ,data = data, subset = train,
                                mtry = ncol(data)-1, importance = TRUE, replace = TRUE, ntree = (i * 10))
  yhat_ite <- predict(bagg_model_cv, newdata = data[-train,])
  rmse_cv_bag[i] = sqrt(mean((yhat_ite - data$Admit[-train])^2))
}


# pick the best number of trees
nTrees = c(1:10)*10
plot(nTrees, rmse_cv_bag, ylab="test RMSE", xlab="n Trees", main="Choosing the best number of trees", type = 'b')

for(i in 1:10){
  if(rmse_cv_bag[i] == min(rmse_cv_bag)){
    abline(v = i * 10, col = "blue", lty = "dashed")
  }
}

nTrees_min_bagg <- 0
for(i in 1:10){
  if(rmse_cv_bag[i] == min(rmse_cv_bag)){
    nTrees_min_bagg <- i * 10
  }
}


# validation
bagg_model <- randomForest(Admit ~ . ,data = data, subset = train,
                              mtry = ncol(data)-1, importance = TRUE, replace = TRUE, ntree = nTrees_min_bagg)
yhat <- predict(bagg_model, newdata = data[-train,])
res_bagg = yhat - data$Admit[-train]
test_RMSE_bagging <- sqrt(mean((res_bagg)^2))
importance(bagg_model)
bagg_model



# residuals
bptest(bagg_model) 
shapiro.test(res_bagg) 
par(mfrow = c(1,4))
hist(res_bagg,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(res_bagg, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(res_bagg),digits = 4),
                                      "- var:", round(var(res_bagg),digits = 4)))
abline(c(0,0),c(0,length(res_bagg)), col= "red", lwd = 2)

boxplot(res_bagg, ylab = "Residuals", main = "Outliers")$out

qqnorm(res_bagg, main='Residuals')
qqline(res_bagg)


################################################################################
#########################################  Random forest
################################################################################

library(randomForest)

set.seed(1)

# Cross validation n trees and selection predictios
forest_cv_models <- list()
for(x in 1:6){
  forest_model_ite <- randomForest(Admit ~ . ,data = data, subset = train,
                                   mtry = x, importance = TRUE, replace = TRUE, ntree = 150)
  forest_cv_models <- append(forest_cv_models, list(forest_model_ite))
}


# create matrix with nReg per split, RMSE min, nTrees
rmse_min_reg <- matrix(list(), nrow=6, ncol=3)
for(x in 1:6){
  for(i in 1:150){
    if(forest_cv_models[[x]]$mse[i] == min(forest_cv_models[[x]]$mse)){
      rmse_min_reg[[x, 1]] <- round(sqrt(min(forest_cv_models[[x]]$mse)), digits = 4)
      rmse_min_reg[[x, 2]] <- i
      
      
      
      forest_model_ite <- randomForest(Admit ~ . ,data = data, subset = train,
                                   mtry = x, importance = TRUE, ntree = i)
      
      res_rdn = yhat - data$Admit[-train]
      yhat <- predict(forest_model_ite, newdata = data[-train,])
      test_RMSE_forest <- sqrt(mean(res_rdn^2))
      
      
      
      rmse_min_reg[[x, 3]] <- round(test_RMSE_forest, digits = 4)
    }
  }
}






par(mfrow = c(3,2))
plot(sqrt(forest_cv_models[[1]]$mse),type = 'b', ylab = "RMSE", xlab = "number of trees", main = "Random forest with 1 regressor per split")
abline(v = index_min[1], col = "blue", lty = "dashed")
plot(sqrt(forest_cv_models[[2]]$mse),type = 'b', ylab = "RMSE", xlab = "number of trees", main = "Random forest with 2 regressor per split")
abline(v = index_min[2], col = "blue", lty = "dashed")
plot(sqrt(forest_cv_models[[3]]$mse),type = 'b', ylab = "RMSE", xlab = "number of trees", main = "Random forest with 3 regressor per split")
abline(v = index_min[3], col = "blue", lty = "dashed")
plot(sqrt(forest_cv_models[[4]]$mse),type = 'b', ylab = "RMSE", xlab = "number of trees", main = "Random forest with 4 regressor per split")
abline(v = index_min[4], col = "blue", lty = "dashed")
plot(sqrt(forest_cv_models[[5]]$mse),type = 'b', ylab = "RMSE", xlab = "number of trees", main = "Random forest with 5 regressor per split")
abline(v = index_min[5], col = "blue", lty = "dashed")
plot(sqrt(forest_cv_models[[6]]$mse),type = 'b', ylab = "RMSE", xlab = "number of trees", main = "Random forest with 6 regressor per split")
abline(v = index_min[6], col = "blue", lty = "dashed")



# best model validation
forest_model <- randomForest(Admit ~ . ,data = data, subset = train,
                             mtry = 2, importance = TRUE, ntree = 78)

res_rdn = yhat - data$Admit[-train]
yhat <- predict(forest_model, newdata = data[-train,])
test_RMSE_forest <- round(sqrt(mean(res_rdn^2)), digits = 4) 

forest_model
importance(forest_model)




# residuals
bptest(forest_model) 
shapiro.test(res_rdn)
par(mfrow = c(1,4))
hist(res_rdn,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(res_rdn, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(res_rdn),digits = 4),
                                      "- var:", round(var(res_rdn),digits = 4)))
abline(c(0,0),c(0,length(res_rdn)), col= "red", lwd = 2)

boxplot(res_rdn, ylab = "Residuals", main = "Outliers")$out

qqnorm(res_rdn, main='Residuals')
qqline(res_rdn)


################################################################################
#########################################  Boosting
################################################################################

install.packages("gbm")
library(gbm)
set.seed(1)



ntree = 6000; 
rmse_cv_boost <- rep(0,3)
boosting_cv_models <- list()
for(i in 1:3){
  boost_model_ite <- gbm(Admit ~ . , data = data[train,], 
                         distribution = "gaussian" , n.trees = ntree,
                         interaction.depth = i, shrinkage = 0.001, verbose = F)
  yhat_ite <- predict(boost_model_ite, newdata = data[-train,], n.trees = ntree)
  rmse_cv_boost[i] <- round(sqrt(mean((yhat_ite - data$Admit[-train])^2)), digits = 5)
  boosting_cv_models <- append(boosting_cv_models, list(boost_model_ite))
}



# plotting the 3 models using the same learning rate but different mumber of regressors
# per split
par(mfrow = c(1,1))
plot(sqrt(boosting_cv_models[[1]]$train.error),  ylab = "train RMSE", xlab = "number of trees", main = "Boosting model with different deapths", col = "#5EC284", lwd=1)
par(new=TRUE)
plot(sqrt(boosting_cv_models[[2]]$train.error), axes=FALSE,  ann=FALSE, col = "#738DF9")
par(new=TRUE)
plot(sqrt(boosting_cv_models[[3]]$train.error),  axes=FALSE,  ann=FALSE, col = "#D28124")
legend(5000, 0.12, legend=c("depth = 1", "depth = 2", "depth = 3"), fill =c("#5EC284", "#738DF9", "#D28124"))


# min depth
min_depth <- rep(0,1)
for(i in 1:3){
  if(rmse_cv_boost[i] == min(rmse_cv_boost)){
    min_depth <- i
  }
}


# validation
boost_model <- gbm(Admit ~ . , data = data[train,], 
                   distribution = "gaussian" , n.trees = ntree,
                   interaction.depth = min_depth, shrinkage = 0.001 , verbose = F)
yhat <- predict(boost_model, newdata = data[-train,], n.trees = ntree)
res_boostin = yhat - data$Admit[-train]
test_RMSE_boost <- round(sqrt(mean((res_boostin)^2)), digits = 5)


# residuals
shapiro.test(res_boostin)
par(mfrow = c(1,4))
hist(res_boostin,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(res_boostin, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(res_boostin),digits = 4),
                                      "- var:", round(var(res_boostin),digits = 4)))
abline(c(0,0),c(0,length(res_boostin)), col= "red", lwd = 2)

boxplot(res_boostin, ylab = "Residuals", main = "Outliers")$out

qqnorm(res_boostin, main='Residuals')
qqline(res_boostin)


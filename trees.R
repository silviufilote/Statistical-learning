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


# Outliers detaction
detect_outlier <- function(x) {
  Quantile1 <- quantile(x, probs=.25)
  Quantile3 <- quantile(x, probs=.75)
  IQR = Quantile3 - Quantile1
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}


remove_outlier <- function(dataframe,columns=names(dataframe)) {
  for (col in columns) {
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  print("Remove outliers")
  print(dataframe)
}


# Load dataset and setup
data <- read.csv("C:\\Users\\fsilv\\Dropbox\\UNI\\stistics\\1.csv")
data <- data.frame(data[,2:9])
colnames(data)[1] = "GRE"
colnames(data)[2] = "TOEFL"
colnames(data)[3] = "UniRatings"
colnames(data)[8] = "Admit"

data <- remove_outlier(data, c('GRE', 'TOEFL', 'UniRatings', 'SOP', "LOR", 'CGPA', "Research"))
data <- data %>% relocate(CGPA, .after=TOEFL)
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

tree_model <- tree( Admit ~ . , data, split = "gini")
summary(tree_model)
plot(tree_model)
text(tree_model,pretty = 0)


model_tree_cv <- cv.tree(tree_model,  FUN = prune.tree)

# dev = cross-validation error rate
# size = number of terminal nodes of each tree -> complexity
plot(model_tree_cv$size, model_tree_cv$dev, ylab = "cross-validation error rate", xlab = "Complexity = foglie")

best = min(model_tree_cv$size[model_tree_cv$dev == min(model_tree_cv$dev)])
prune <- prune.tree(tree_model, best = best)
summary(prune)
plot(prune)
text(prune, pretty = 0)



# validation and MSE
tree_model_train <- tree( Admit ~ . , data, subset = train, split = "gini")
train_prune <- prune.tree(tree_model_train, best = best)
summary(tree_model_train)

pred_value <- predict(tree_model_train, newdata = data[-train,])
# table(pred_value, data$Admit[-train])


res = data$Admit[-train] - pred_value
test_MSE_tree = mean((res)^2)
bptest(tree_model_train) # non sono omoschedastici
shapiro.test(res) # non sono normali


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
grid()

# utilizzo 500 alberi
bagg_test$ntree

# validation
yhat <- predict(bagg_test, newdata = data[-train,])
test_MSE_bagging <- mean((yhat - data$Admit[-train])^2)

bptest(bagg_model) # non sono omoschedastici
shapiro.test((yhat - data$Admit[-train])) # non sono normali 
importance(bagg_test)


# fare crossvalidazione sul numero di alberi
mse_cv_bag = rep(0,100)
for (i in 1:100){
  bagg_model_cv <- randomForest(Admit ~ . ,data = data, subset = train,
                                mtry = ncol(data)-1, importance = TRUE, replace = TRUE, ntree = (i * 10))
  yhat_ite <- predict(bagg_model_cv, newdata = data[-train,])
  mse_cv_bag[i] = mean((yhat_ite - data$Admit[-train])^2)
}

treesX = c(1:100)*10
plot(treesX, mse_cv_bag, type = 'b')
grid()


nTrees_bagg <- 0
for(i in 1:100){
  if(mse_cv_bag[i] == min(mse_cv_bag)){
    nTrees_bagg <- i * 10
  }
}

################################################################################
#########################################  Random forest
################################################################################


forest_model <- randomForest(Admit ~ . ,data = data, subset = train,
                             mtry = floor(sqrt(ncol(data)-1)), importance = TRUE)


plot(forest_model,type = 'b',col="green",pch = "+")
par(new=TRUE)
plot(bagg_test,type = 'b',col="red",pch='o')

yhat <- predict(forest_model, newdata = data[-train,])
test_MSE_forest <- mean((yhat - data$Admit[-train])^2)

forest_model
importance(forest_model)
bptest(forest_model) 
shapiro.test((yhat - data$Admit[-train])) # non sono normali 

# Cross validation n trees and selection predictios

mse_cv_forest <- matrix(nrow=6, ncol=100)
for(x in 1:6){
  for (i in 1:100){
    forest_model_ite <- randomForest(Admit ~ . ,data = data, subset = train,
                                    mtry = x, importance = TRUE, replace = TRUE, ntree = (i * 10))
    
    forest_cv_models[[x,i]] <- forest_model_ite
    yhat_ite <- predict(forest_model_ite, newdata = data[-train,])
    mse_cv_forest[[x,i]] = mean((yhat_ite - data$Admit[-train])^2)
  }
}

par(mfrow = c(3,2))
plot(treesX, mse_cv_forest[1,], type = 'b')
plot(treesX, mse_cv_forest[2,], type = 'b')
plot(treesX, mse_cv_forest[3,], type = 'b')
plot(treesX, mse_cv_forest[4,], type = 'b')
plot(treesX, mse_cv_forest[5,], type = 'b')
plot(treesX, mse_cv_forest[6,], type = 'b')


nTrees_forest <- matrix(nrow=2, ncol=6)

for(x in 1:6){
  for(i in 1:100){
    if(mse_cv_forest[i] == min(mse_cv_forest[x,])){
      nTrees_forest[1,x] <- (i * 10)
      nTrees_forest[2,x] <- mse_cv_forest[i]
    }
  }
}


test_MSE_forest = min(mse_cv_forest)

################################################################################
#########################################  Random forest
################################################################################

install.packages("gbm")
library(gbm)
set.seed(1)

ntree = 5000; 
boost_model <- gbm(Admit ~ . , data = data[train,], 
                   distribution = "gaussian" , n.trees = ntree,
                   interaction.depth = 4, shrinkage = 0.01 , verbose = F)
boost_model

yhat <- predict(boost_model, newdata = data[-train,], n.trees = ntree)

# test mse migliore dell'altro
test_MSE_boost <- mean((yhat - data$Admit[-train])^2)


mse_cv_boost <- matrix(nrow=3, ncol=40)
leaning_rate = c(0.01, 0.1, 1)
span_x = 1:40


for(x in 1:3){
  for(i in span_x){
    ntree = (1000*i)
    boost_model_ite <- gbm(Admit ~ . , data = data[train,], 
                       distribution = "gaussian" , n.trees = ntree,
                       interaction.depth = 4, shrinkage = leaning_rate[x] , verbose = F)
    yhat_ite <- predict(boost_model_ite, newdata = data[-train,], n.trees = ntree)
    mse_cv_boost[[x,i]] <- mean((yhat_ite - data$Admit[-train])^2)
  }
}


par(mfrow = c(1,3))
plot(span_x*1000, mse_cv_boost[1,], xlab = "trees", ylab = "test mse", type = 'b')
plot(span_x*1000, mse_cv_boost[2,], xlab = "trees", ylab = "test mse", type = 'b')
plot(span_x*1000, mse_cv_boost[3,], xlab = "trees", ylab = "test mse", type = 'b')


nTrees_forest <- matrix(nrow=2, ncol=6)

for(x in 1:6){
  for(i in 1:100){
    if(mse_cv_forest[i] == min(mse_cv_forest[x,])){
      nTrees_forest[1,x] <- (i * 10)
      nTrees_forest[2,x] <- mse_cv_forest[i]
    }
  }
}


test_MSE_boost = min(mse_cv_boost)



################################################################################
#########################################  Random forest
################################################################################


# fare crossvalidazione sul numero di tree -> scegliere il migliore
# fare crossvalidazione sul numero di regressori da considerare ad ogni split




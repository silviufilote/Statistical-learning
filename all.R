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

# ctrl + L => to clear console

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
# library(olsrr)

################################################################################
######################################### Dataset 
################################################################################

data <- read.csv("C:\\Users\\fsilv\\Dropbox\\UNI\\stistics\\1.csv")
data <- data.frame(data[,2:9])
colnames(data)[1] = "GRE"
colnames(data)[2] = "TOEFL"
colnames(data)[3] = "UniRatings"
colnames(data)[8] = "Admit"
summary(data)
head(data)


# no nan colums
colSums(is.na(data))



################################################################################
#########################################  Istogrammi
################################################################################

par(mfrow = c(3,3))

hist(data$GRE, 40,
     xlab = "GRE",
     main = "GRE SCORE distribution") 


hist(data$TOEFL, 40,
     xlab = "TOEFL",
     main = "TOEFL SCORE distribution ") 


hist(data$UniRatings, 40,
     xlab = "University Rating",
     main = "University Rating distribution")


hist(data$SOP, 40,
     xlab = "SOP",
     main = "SOP distribution")


hist(data$LOR, 40,
     xlab = "LOR",
     main = "LOR distribution")

hist(data$CGPA, 40,
     xlab = "CGPA",
     main = "CGPA distribution")


hist(data$Research, 40,
     xlab = "Research",
     main = "Research distribution")


hist(data$Admit, 40,
     xlab = "Chance of Admit",
     main = "Chance of Admit distribution")




################################################################################
#########################################  Analisi univariata
################################################################################


par(mfrow = c(3,3))
plot(data$GRE, 
     data$Admit, 
     main = "GRE score vs Chance of Admit", 
     xlab = "GRE score", 
     ylab = "Chance of Admit")
abline(lm(data$Admit ~ data$GRE, data = mtcars), col = "blue")

plot(data$TOEFL, 
     data$Admit,  
     main = "TOEFL score vs Chance of Admit", 
     xlab = "TOEFL score", 
     ylab = "Chance of Admit")
abline(lm(data$Admit ~ data$TOEFL, data = mtcars), col = "blue")

plot(data$UniRatings, 
     data$Admit,  
     main = "UniRatings vs Chance of Admit", 
     xlab = "UniRatings", 
     ylab = "Chance of Admit")
abline(lm(data$Admit ~ data$UniRatings, data = mtcars), col = "blue")

plot(data$SOP, 
     data$Admit,  
     main = "SOP vs Chance of Admit", 
     xlab = "SOP", 
     ylab = "Chance of Admit")
abline(lm(data$Admit ~ data$SOP, data = mtcars), col = "blue")

plot(data$LOR, 
     data$Admit,  
     main = "LOR vs Chance of Admit", 
     xlab = "LOR", 
     ylab = "Chance of Admit")
abline(lm(data$Admit ~ data$LOR, data = mtcars), col = "blue")

plot(data$CGPA, 
     data$Admit,  
     main = "CGPA vs Chance of Admit", 
     xlab = "CGPA", 
     ylab = "Chance of Admit")
abline(lm(data$Admit ~ data$CGPA, data = mtcars), col = "blue")

plot(data$Research, 
     data$Admit, 
     main = "Research vs Chance of Admit", 
     xlab = "Research", 
     ylab = "Chance of Admit")
abline(lm(data$Admit ~ data$Research, data = mtcars), col = "blue")



################################################################################
#########################################  pairs-plot
################################################################################

library(ggplot2)
library(GGally)

par(mfrow = c(1,1))
cor_scores <- cor(data)
corrplot(cor_scores, method="number", type = 'lower')

ggpairs(data)

################################################################################
############################### Outliers
################################################################################

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

data <- remove_outlier(data, c('GRE', 'TOEFL', 'UniRatings', 'SOP', "LOR", 'CGPA', "Research"))


# crate multiple dataset to fit later
data <- data %>% relocate(CGPA, .after=TOEFL)


# after standirzation lets check the outliers
par(mfrow = c(1,1))
boxplot(scale(data))$out



################################################################################
# fit linear model: using whole data ===> no dummy
################################################################################


# fit the model
datanf <- data
data <- data %>% 
  mutate_at(vars(Research, LOR, SOP, UniRatings), as.factor)

lm_fit <- lm(Admit ~ ., data = data) 
lm_fit <- update(lm_fit, ~ . - SOP - UniRatings - LOR)
summary(lm_fit)




############### Prove
shapiro.test(lm_fit$residuals)
anova(lm_fit)
confint(lm_fit, level = 0.98)



################################################################################
# fit linear model: using dummy
################################################################################


lm_fit_dummy_all <- lm(Admit ~ ., data = data) 
summary(lm_fit_dummy_all)



# Step-wise Regression -> best model based on AIC value.

lm_fit_dummy_none <- lm(formula = Admit ~ 1, data = data)
model_stepwise <- step(
  object = lm_fit_dummy_none,
  direction = "both",
  scope = list(upper = lm_fit_dummy_all),
  trace = FALSE)
summary(model_stepwise)


# model performance
performance <- compare_performance(lm_fit, lm_fit_dummy_all, model_stepwise)
as.data.frame(performance)



################################################################################
#########################################  Tests :
######################################### - normality of residuals
######################################### - homoschedasticity of residuals
################################################################################


par(mfrow = c(1,1))
boxplot(lm_fit$residuals)$out
var(lm_fit$residuals)


par(mfrow = c(2,2))

hist(lm_fit$residuals,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(lm_fit$residuals, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(lm_fit$residuals),digits = 4),
                                      "- var:", round(var(lm_fit$residuals),digits = 2)))
abline(c(0,0),c(0,length(lm_fit$residuals)), col= "red", lwd = 2)

boxplot(lm_fit$residuals)$out

qqnorm(lm_fit$residuals, main='Residuals')
qqline(lm_fit$residuals)



set.seed(1)

# Shapiro-Wilk hypothesis test:
# H0: Variable is normally distributed
# H1: Variable is not normally distributed
# p-value significative => residuals are not normally distributed
shapiro.test(lm_fit$residuals)
shapiro.test(lm_fit_dummy_all$residuals)
shapiro.test(model_stepwise$residuals)


# Breusch-Pagan hypothesis test:
# H0: Homoscedasticity is present (the residuals are distributed with equal variance)
# H1 : Heteroscedasticity is present (the residuals are not distributed with equal variance)
# p-value significative => residuals arent homoscedastical
bptest(lm_fit)
bptest(lm_fit_dummy_all)
bptest(model_stepwise)



################################################################################
#########################################  Testing betas => Std. Error
################################################################################

library(boot)

set.seed(1)

get_model <- function(data, index){
  boot_train <- data[index,]
  boot_model <- glm(Admit ~ ., data = boot_train)  
  boot_model <- update(boot_model, ~ . - SOP - UniRatings - LOR)
  return (boot_model$coefficients)
}

boot_model <- boot(data = data, statistic = get_model, R = 1000);
boot_model$t0
boot_model




################################################################################
#########################################  Model Improvement
################################################################################


# Transform Target-Var (Arcsin)
admission_transform_y  <- data %>%
  select(Admit, GRE, TOEFL, LOR, CGPA, Research) %>%
  mutate(Admit = asin(sqrt(Admit)))
fit_arcisn <- lm(Admit ~ + GRE + TOEFL + CGPA + Research, admission_transform_y)
summary(fit_arcisn)

shapiro.test(fit_arcisn$residuals) # no normality
bptest(fit_arcisn) #  no homoscehdasticity



# Log10 Regressors
fit_log10 <- lm(asin(sqrt(data$Admit)) ~ log10(GRE) + log10(TOEFL) + log10(CGPA) + Research, data)
summary(fit_log10)
shapiro.test(fit_log10$residuals) # no normality
bptest(fit_log10) #  homoscehdasticity


# Standardized regressors
fit_scale <- lm(asin(sqrt(data$Admit)) ~ scale(GRE) + scale(TOEFL) + scale(CGPA) + Research, data)
summary(fit_scale)
shapiro.test(fit_scale$residuals) # no normality
bptest(fit_scale) #  no homoscehdasticity


# Sqrt regressors
fit_sqrt <- lm(asin(sqrt(data$Admit)) ~ sqrt(GRE) + sqrt(TOEFL) + sqrt(CGPA) + Research, data)
summary(fit_sqrt)
shapiro.test(fit_sqrt$residuals) # no normality
bptest(fit_sqrt) #  homoscehdasticity


performance_model <- compare_performance(lm_fit, fit_arcisn, fit_log10, fit_scale, fit_sqrt)
as.data.frame(performance_model)


################################################################################
#########################################  Validation
################################################################################

library(boot)

summary(lm_fit)
# mean(lm_fit$residuals^2)
total_MSE = anova(lm_fit)['Residuals', 'Mean Sq']


glm_fit <- glm(asin(sqrt(data$Admit)) ~ sqrt(data$GRE) + sqrt(data$TOEFL) + sqrt(data$CGPA) + data$Research, data = data)
model_cv <- cv.glm(data, glm_fit)
LOOCV_MSE <- model_cv$delta
LOOCV_MSE_train <- model_cv$delta[1]
LOOCV_MSE_test <- model_cv$delta[2]


cv.mse = rep(0,3)
res_norm = rep(0,3)
res_homo = rep(0,3)
for (i in 1:3){
  model_ite <- glm(asin(sqrt(data$Admit)) ~  poly(sqrt(data$GRE), i) + poly(sqrt(data$TOEFL), i) + poly(sqrt(data$CGPA), i) + data$Research, data = data)
  cv.mse[i] = cv.glm(data, model_ite)$delta[2]
  results <- shapiro.test(model_ite$residuals)
  res_norm[i] <- results$p.value
  results <- bptest(model_ite)
  res_homo[i] = results$p.value
}

cv.mse
res_norm
res_homo

par(mfrow = c(2,2))
plot(1:3,cv.mse ,type = "b",col = "blue",
     ylab = "TEST MSE",
     xlab = "Flexibility (poly degree)",
     main = "Test error estimation")
plot(1:3,res_norm ,type = "b",col = "blue",
     ylab = "p-value normality res",
     xlab = "Flexibility (poly degree)",
     main = "Normality residuals")
plot(1:3,res_homo ,type = "b",col = "blue",
     ylab = "p-value homoshedasticity res",
     xlab = "Flexibility (poly degree)",
     main = "Homoschedasticity residuals")




################################################################################
#########################################  LASSO 
################################################################################

library(boot)
library(glmnet)

df <- data

for(x in 1:1:(dim(df)[2] - 5)){
  df[paste("l",colnames(df)[x], sep = "")] <- log(df[x])
  df[paste("s",colnames(df)[x], sep = "")] <- (df[x])^2
  df[paste("c",colnames(df)[x], sep = "")] <- (df[x])^3
  df[paste("r",colnames(df)[x], sep = "")] <- sqrt(df[x])
  df[paste("e",colnames(df)[x], sep = "")] <- exp(df[x])
}

x <- model.matrix(Admit ~ . ,df) [, -1]
y <- df$Admit
lambda <- 10^seq(-3, 0,length = 300)

set.seed(1)
train <- sample(nrow(x), floor(nrow(x)*0.5), replace = FALSE)
cv_lasso <- cv.glmnet(x[train, ], y[train], alpha = 1, lambda = lambda );
plot(cv_lasso)


coef(cv_lasso, cv_lasso$lambda.1se)
coef(cv_lasso, cv_lasso$lambda.min)

opt_lambda <- cv_lasso$lambda.min

lasso_model <- glmnet(x[train,], y[train],alpha = 1, lambda = opt_lambda)
fitt_value <- predict(lasso_model, s = opt_lambda, newx=x[-train,])
test_MSE_lasso = mean((y[-train] - fitt_value)^2)

res_lasso = y[-train] - fitt_value



# Testing residuals
shapiro.test(res_lasso)

par(mfrow = c(2,2))

hist(res_lasso,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(res_lasso, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(res_lasso),digits = 4),
                                      "- var:", round(var(res_lasso),digits = 2)))
abline(c(0,0),c(0,length(res_lasso)), col= "red", lwd = 2)

boxplot(res_lasso)$out

qqnorm(res_lasso, main='Residuals')
qqline(res_lasso)


################################################################################
#########################################  Step functions
################################################################################

table(cut(datanf$LOR, 10))

train = sample(dim(datanf)[1], dim(datanf)[1]*0.7, replace = FALSE)

step_fun <- lm(Admit ~ GRE + TOEFL + cut(LOR, 10) + CGPA + cut(Research,2), data = datanf);
summary(step_fun)

step_fun_train <- lm(Admit ~ GRE + TOEFL + cut(LOR, 10) + CGPA + cut(Research,2), data = datanf, subset = train);
summary(step_fun_train)
res_step = (datanf$Admit - predict(step_fun_train, datanf))^2
train_err_SF = mean(res_step[train])
test_MSE_Step = mean(res_step[-train])



# Testing residuals
shapiro.test(step_fun$residuals)
bptest(step_fun_train) 


################################################################################
#########################################  Tree
################################################################################

install.packages("tree")
library(tree)

tree_model <- tree( Admit ~ . , data, split = "gini")
summary(tree_model)
plot(tree_model)
text(tree_model,pretty = 0)


tree_model_train <- tree( Admit ~ . , data, subset = train)
summary(tree_model)

pred_value <- predict(tree_model_train, newdata = data[-train,])
table(pred_value, data$Admit[-train])

test_MSE_tree = mean((data$Admit[-train] - pred_value)^2)
bptest(tree_model_train) 

################################################################################
#########################################  Bagging model
################################################################################
install.packages("randomForest")
library(randomForest)
set.seed(1)

bagg_model <- randomForest(Admit ~ . ,data = data, subset = train,
                           mtry = ncol(data)-1, importance = TRUE, replace = TRUE)
bagg_model = bagg_model$mse
plot(bagg_model)

yhat <- predict(bagg_model, newdata = data[-train,])
plot(yhat,data$Admit[-train])
abline(0,1)

test_MSE_bagging <- mean((yhat - data$Admit[-train])^2)

# fare crossvalidazione sul numero di alberi
bagg_model <- randomForest(Admit ~ . ,data = data, subset = train,
                           mtry = ncol(data)-1, importance = TRUE, replace = TRUE, ntree = 100)
bagg_model
importance(bagg_model)


bptest(bagg_model) 

################################################################################
#########################################  Random forest
################################################################################


forest_model <- randomForest(Admit ~ . ,data = data, subset = train,
                             mtry = floor(sqrt(ncol(data)-1)), importance = TRUE, ntree = 100)

forest_model

plot(forest_model,type = 'b',col="green",pch = "+")
par(new=TRUE)
plot(bagg_model,type = 'b',col="red",pch='o')

yhat <- predict(forest_model, newdata = data[-train,])
test_MSE_rand <- mean((yhat - data$Admit[-train])^2)
importance(forest_model)

bptest(forest_model) 

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


################################################################################
#########################################  Random forest
################################################################################


# fare crossvalidazione sul numero di tree -> scegliere il migliore
# fare crossvalidazione sul numero di regressori da considerare ad ogni split











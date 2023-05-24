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
library(boot)
library(glmnet)


# dataset cleaning 
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

data <- data %>% relocate(CGPA, .after=TOEFL)


################################################################################
############################### Train model
################################################################################

set.seed(1)
train = sample(dim(data)[1], dim(data)[1]*0.7, replace = FALSE)

df = data[train,]

data <- data %>% 
  mutate_at(vars(Research, LOR, SOP, UniRatings), as.factor)


################################################################################
############################### fit linear model: using whole data
################################################################################


lm_model <- lm(Admit ~ ., data = data, subset = train) 
summary(lm_model,digits=3)


lm_model <- update(lm_model, ~ . - SOP - UniRatings)
summary(lm_model,digits=3)


#shapiro.test(lm_fit$residuals)
# anova(lm_fit_train)
## confint(lm_fit, level = 0.98)
## total_MSE = anova(lm_fit)['Residuals', 'Mean Sq']


# analisys res -> train set
par(mfrow = c(1,4))
hist(lm_model$residuals,40,
     xlab = "Residual",
     main = "histogram residuals") 

plot(lm_model$residuals, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(lm_model$residuals),digits = 4),
                                      "- var:", round(var(lm_model$residuals),digits = 4)))
abline(c(0,0),c(0,length(lm_model$residuals)), col= "red", lwd = 2)

boxplot(lm_model$residuals, ylab = "Residuals", main = "Outliers")$out

qqnorm(lm_model$residuals, main='Residuals')
qqline(lm_model$residuals)

shapiro.test(lm_model$residuals) # non normality
bptest(lm_model) # non homoscehdasticity



# Validation 
yhat <- predict(lm_model, newdata = data[-train,])
lm_RMSE <- round(sqrt(mean((yhat - data$Admit[-train])^2)), digits = 4)


for(i in 1:dim(data.frame(yhat))[1]){
  if(yhat[i] < 0 || yhat[i] > 1){
    print(yhat[i])
  }
}



# Cross-validation
model_cv = glm(Admit ~ .  - SOP - UniRatings, data = data)
model_cv = cv.glm(data, model_cv)
lm_rMSE_train_LOOCV <- round(sqrt(model_cv$delta[1]), digits = 4)
lm_rMSE_test_LOOCV <- round(sqrt(model_cv$delta[2]), digits = 4)


################################################################################
###################################### Testing betas lm_fit model => Std. Error
################################################################################


library(boot)

set.seed(1)

# train <- sample(dim(df)[1], dim(df)[1]*0.7, replace = FALSE)
# df = data[train,]

df <- df %>% 
  mutate_at(vars(Research, LOR, SOP, UniRatings), as.factor)


get_model <- function(dataset, index){
  boot_data <- dataset[index,]
  boot_model <- glm(Admit ~ GRE + TOEFL + CGPA + LOR + Research, data = boot_data)  
  return (boot_model$coefficients)
}

boot_model <- boot(data = data, statistic = get_model, R = 1000) 
boot_model$t0
boot_model

################################################################################
############################### Correlazione covariate residui 
################################################################################


par(mfrow = c(2,2))
plot(lm_fit_train$model$GRE, 
     lm_fit_train$residuals, 
     main = "GRE score vs Res", 
     xlab = "GRE score", 
     ylab = "Res")
abline(lm(lm_fit_train$residuals ~ lm_fit_train$model$GRE, data = data), col = "blue")

plot(lm_fit_train$model$TOEFL, 
     lm_fit_train$residuals,  
     main = "TOEFL score vs Res", 
     xlab = "TOEFL score", 
     ylab = "Res")
abline(lm(lm_fit_train$residuals ~ lm_fit_train$model$TOEFL, data = data), col = "blue")

plot(lm_fit_train$model$CGPA, 
     lm_fit_train$residuals,  
     main = "CGPA vs Res", 
     xlab = "CGPA", 
     ylab = "Res")
abline(lm(lm_fit_train$residuals ~ lm_fit_train$model$CGPA, data = data), col = "blue")

plot(lm_fit_train$model$Research, 
     lm_fit_train$residuals,  
     main = "Research vs Res", 
     xlab = "Research", 
     ylab = "Res")
abline(lm(lm_fit_train$residuals ~ lm_fit_train$model$Research, data = data), col = "blue")


################################################################################
############################### Stepwise regression method
################################################################################

set.seed(1)

lm_fit_all <- lm(Admit ~ . + exp(GRE) + exp(TOEFL) + exp(CGPA) + log10(GRE) + log10(TOEFL) + log10(CGPA) + poly(GRE,3) + poly(TOEFL,3) + poly(CGPA,3) + sqrt(GRE) + sqrt(TOEFL) + sqrt(CGPA), data = data, subset = train) 
summary(lm_fit_all)


# Step-wise Regression -> best model based on AIC value.
lm_fit_none <- lm(formula = Admit ~ 1, data = data, subset = train)

stepwise_model <- step(
  object = lm_fit_all,
  direction = "backward",
  trace = FALSE)
summary(stepwise_model)


# model performance
performance <- compare_performance(lm_fit_train, lm_fit_all, stepwise_model)
as.data.frame(performance)




 # analisys res -> all dataset
shapiro.test(stepwise_model$residuals)
bptest(stepwise_model)

par(mfrow = c(1,4))
hist(stepwise_model$residuals,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(stepwise_model$residuals, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(stepwise_model$residuals),digits = 4),
                                      "- var:", round(var(stepwise_model$residuals),digits = 4)))
abline(c(0,0),c(0,length(stepwise_model$residuals)), col= "red", lwd = 2)

boxplot(stepwise_model$residuals, ylab = "Residuals", main = "Outliers")$out

qqnorm(stepwise_model$residuals, main='Residuals')
qqline(stepwise_model$residuals)



# Validation 
yhat <- predict(stepwise_model, newdata = data[-train,])
step_RMSE_test_val <- round(sqrt(mean((yhat - data$Admit[-train])^2)), digits = 4)



for(i in 1:dim(data.frame(yhat))[1]){
  if(yhat[i] < 0 || yhat[i] > 1){
    print(yhat[i])
  }
}


# Cross-validation
model_cv <- glm(Admit ~ log10(GRE) + log10(TOEFL) + log10(CGPA) + exp(CGPA) + poly(CGPA, 3) + LOR + Research + UniRatings, data = data)
model_cv <- cv.glm(data, model_cv)
step_MSE_train_LOOCV <- round(sqrt(model_cv$delta[1]), digits = 4)
step_MSE_test_LOOCV <- round(sqrt(model_cv$delta[2]), digits = 4)




################################################################################
################################### Testing betas stepwise model => Std. Error
################################################################################

library(boot)

set.seed(1)

get_model <- function(data, index){
  boot_train <- data[index,]
  boot_model <- glm(Admit ~ CGPA + GRE + LOR + Research + TOEFL + UniRatings, data = boot_train)  
  boot_model$coefficients
}

boot_result  <- boot(data = data, statistic = get_model, R = 1000);
boot_result$t0
boot_result 

# Calculate the bootstrap mean
bootstrap_mean <- mean(boot_result$t)

# Obtain bootstrap confidence intervals
bootstrap_ci <- boot.ci(boot_result, type = "basic")


################################################################################
#########################################  LASSO 
################################################################################

set.seed(1)
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

df[paste("e",colnames(df)[x], sep = "")] <- exp(df[x])


x <- model.matrix(Admit ~ . ,df) [, -1]
y <- df$Admit
lambda <- 10^seq(-3, 0,length = 300)

set.seed(1)
cv_lasso <- cv.glmnet(x[train, ], y[train], alpha = 1, lambda = lambda );

par(mfrow = c(1,1))
plot(cv_lasso)


coef(cv_lasso, cv_lasso$lambda.1se)
coef(cv_lasso, cv_lasso$lambda.min)


opt_lambda <- cv_lasso$lambda.1se



lasso_model <- glmnet(x[train,], y[train], alpha = 1, lambda = opt_lambda)
fitt_value <- predict(lasso_model, s = opt_lambda, newx=x[-train,])
lasso_RMSE_test = round(sqrt(mean((y[-train] - fitt_value)^2)), digits = 4)

res_lasso = y[-train] - fitt_value
coef(lasso_model, opt_lambda)


# Testing residuals
shapiro.test(res_lasso) # non normality
#bptest() # vale solo nel caso lineare
lasso_Rsquadre = lasso_model$dev.ratio

par(mfrow = c(1,4))

hist(res_lasso,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(res_lasso, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(res_lasso),digits = 4),
                                      "- var:", round(var(res_lasso),digits = 4)))
abline(c(0,0),c(0,length(res_lasso)), col= "red", lwd = 2)

boxplot(res_lasso, ylab = "Residuals", main = "Outliers", )$out

qqnorm(res_lasso, main='Residuals')
qqline(res_lasso)

model_cv_lasso = glm(Admit ~ (GRE)^2 + (GRE)^3 + TOEFL + (TOEFL)^2 + (TOEFL)^3 + CGPA + Research, data = data)
model_cv_lasso = cv.glm(data, model_cv_lasso)
model_imp_rMSE_lasso_train <- round(sqrt(model_cv_lasso$delta[1]), digits = 4)
model_imp_rMSE_lasso_test <- round(sqrt(model_cv_lasso$delta[2]), digits = 4)




################################################################################
#########################################  Model Improvement
################################################################################


lm_model_imp <- lm( asin(sqrt(Admit)) ~ ., data = data, subset = train) 
lm_model_imp <- update(lm_model_imp, ~ . - SOP - UniRatings)
summary(lm_model_imp)

yhat <- predict(lm_model_imp, newdata = data[-train,])
step_RMSE_test_val <- round(sqrt(mean((sin(yhat^2) - data$Admit[-train])^2)), digits = 4)

shapiro.test(lm_model_imp$residuals) # no normality
bptest(lm_model_imp) 


model_cv_model_imp = glm( asin(sqrt(Admit)) ~ . - SOP - UniRatings, data = data)
model_cv_model_imp = cv.glm(data, model_cv_model_imp)
model_imp_rMSE_train_LOOCV <- round(sqrt(model_cv_model_imp$delta[1]), digits = 4)
model_imp_rMSE_test_LOOCV <- round(sqrt(model_cv_model_imp$delta[2]), digits = 4)



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
bptest(fit_arcisn) #  homoscehdasticity


# Log10 Regressors
fit_log10 <- lm(asin(sqrt(data$Admit)) ~ log10(GRE) + log10(TOEFL) + log10(CGPA) + Research, data)
summary(fit_log10)
shapiro.test(fit_log10$residuals) # no normality
bptest(fit_log10) #  homoscehdasticity


# Standardized regressors
fit_scale <- lm(asin(sqrt(data$Admit)) ~ scale(GRE) + scale(TOEFL) + scale(CGPA) + Research, data)
summary(fit_scale)
shapiro.test(fit_scale$residuals) # no normality
bptest(fit_scale) #  homoscehdasticity


# Sqrt regressors
fit_sqrt <- lm(asin(sqrt(data$Admit)) ~ sqrt(GRE) + sqrt(TOEFL) + sqrt(CGPA) + Research, data)
summary(fit_sqrt)
shapiro.test(fit_sqrt$residuals) # no normality
bptest(fit_sqrt) #  homoscehdasticity


performance_model <- compare_performance(lm_fit, fit_arcisn, fit_log10, fit_scale, fit_sqrt)
as.data.frame(performance_model)






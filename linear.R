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


# clearing environment 
rm(list = ls())
graphics.off()


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


################################################################################
############################### fit linear model: using whole data
################################################################################


data <- data %>% 
  mutate_at(vars(Research, LOR, SOP, UniRatings), as.factor)
lm_fit <- lm(Admit ~ ., data = data) 

lm_fit <- update(lm_fit, ~ . - SOP - UniRatings - LOR)
summary(lm_fit)




# Prove
shapiro.test(lm_fit$residuals)
anova(lm_fit)
confint(lm_fit, level = 0.98)
total_MSE = anova(lm_fit)['Residuals', 'Mean Sq']




# analisys res -> all dataset
shapiro.test(lm_fit$residuals)
bptest(lm_fit)

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




# Validation 
lm_fit_validation <- lm(Admit ~ .  - SOP - UniRatings - LOR, data = data, subset = train)
yhat <- predict(lm_fit_validation, newdata = data[-train,])
lm_MSE_test_val <- mean((yhat - data$Admit[-train])^2)




# Cross-validation
model_cv = glm(Admit ~ .  - SOP - UniRatings - LOR, data = data)
model_cv = cv.glm(data, model_cv)
lm_MSE_train_LOOCV <- model_cv$delta[1]
lm_MSE_test_LOOCV <- model_cv$delta[2]


################################################################################
###################################### Testing betas lm_fit model => Std. Error
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
############################### Stepwise regression method
################################################################################

lm_fit_all <- lm(Admit ~ ., data = data) 
summary(lm_fit_all)




# Step-wise Regression -> best model based on AIC value.
lm_fit_none <- lm(formula = Admit ~ 1, data = data)
model_stepwise <- step(
  object = lm_fit_none,
  direction = "both",
  scope = list(upper = lm_fit_all),
  trace = FALSE)
summary(model_stepwise)




# model performance
performance <- compare_performance(lm_fit, lm_fit_all, model_stepwise)
as.data.frame(performance)




# analisys res -> all dataset
shapiro.test(model_stepwise$residuals)
bptest(model_stepwise)

par(mfrow = c(2,2))
hist(model_stepwise$residuals,40,
     xlab = "Residual",
     main = "Distribuzione empirica dei residui") 

plot(model_stepwise$residuals, pch = "o", col = "blue" ,
     ylab = "Residual", main = paste0("Residual plot - mean:",round(mean(model_stepwise$residuals),digits = 4),
                                      "- var:", round(var(model_stepwise$residuals),digits = 2)))
abline(c(0,0),c(0,length(model_stepwise$residuals)), col= "red", lwd = 2)

boxplot(model_stepwise$residuals)$out

qqnorm(model_stepwise$residuals, main='Residuals')
qqline(model_stepwise$residuals)



# Validation 
lm_fit_none <- lm(formula = Admit ~ 1, data = data, subset = train)
lm_fit_all <- lm(Admit ~ ., data = data, subset = train)
model_stepwise <- step(
  object = lm_fit_none,
  direction = "both",
  scope = list(upper = lm_fit_all),
  trace = FALSE)
yhat <- predict(model_stepwise, newdata = data[-train,])
step_MSE_test_val <- mean((yhat - data$Admit[-train])^2)




# Cross-validation
model_cv <- glm(Admit ~ CGPA + GRE + LOR + Research + TOEFL + UniRatings, data = data)
model_cv <- cv.glm(data, model_cv)
step_MSE_train_LOOCV <- model_cv$delta[1]
step_MSE_test_LOOCV <- model_cv$delta[2]




################################################################################
################################### Testing betas stepwise model => Std. Error
################################################################################

library(boot)

set.seed(1)

get_model <- function(data, index){
  boot_train <- data[index,]
  boot_model <- glm(Admit ~ CGPA + GRE + LOR + Research + TOEFL + UniRatings, data = boot_train)  
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
cv_lasso <- cv.glmnet(x[train, ], y[train], alpha = 1, lambda = lambda );
plot(cv_lasso)


coef(cv_lasso, cv_lasso$lambda.1se)
coef(cv_lasso, cv_lasso$lambda.min)

opt_lambda <- cv_lasso$lambda.min

lasso_model <- glmnet(x[train,], y[train],alpha = 1, lambda = opt_lambda)
fitt_value <- predict(lasso_model, s = opt_lambda, newx=x[-train,])
lasso_MSE_test = mean((y[-train] - fitt_value)^2)

res_lasso = y[-train] - fitt_value


# Testing residuals
shapiro.test(res_lasso) # non normality
bptest(fit_scale) #  homoscehdasticity


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














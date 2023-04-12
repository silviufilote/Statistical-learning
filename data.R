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
     main = "Distribuzione GRE SCORE") 


hist(data$TOEFL, 40,
     xlab = "TOEFL",
     main = "Distribuzione TOEFL SCORE") 


hist(data$UniRatings, 40,
     xlab = "University Rating",
     main = "Distribuzione University Rating")


hist(data$SOP, 40,
     xlab = "SOP",
     main = "Distribuzione SOP")


hist(data$LOR, 40,
     xlab = "LOR",
     main = "Distribuzione LOR")

hist(data$CGPA, 40,
     xlab = "CGPA",
     main = "Distribuzione CGPA")


hist(data$Research, 40,
     xlab = "Research",
     main = "Distribuzione Research")


hist(data$Admit, 40,
     xlab = "Chance of Admit",
     main = "Distribuzione Chance of Admit")




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


remove_outlier <- function(dataframe) {
  columns <- names(dataframe)
  for (col in columns) {
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  print("Remove outliers")
  print(dataframe)
}

# after standirzation lets check the outliers
par(mfrow = c(1,1))
boxplot(scale(data))$out
data <- remove_outlier(data)
data <- remove_outlier(data)


par(mfrow = c(1,1))
boxplot(scale(data))$out



# crate multiple dataset to fit later
data <- data %>% relocate(CGPA, .after=TOEFL)




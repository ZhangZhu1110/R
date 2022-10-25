setwd("C:/Users/zhang/Desktop/R")
library(readxl)
library(ggplot2)
library(skedastic)
library(lmtest)

# The assumption of linear regression:
# Linearity: the relationship between the independent and the dependent variables need to be linear
# No or little multicollinearity: multicollinearity occurs when the independent variables are to high correlated with each other
# No auto-correlation: auto-correlation occurs when the residuals are not independent from each other
# Homoscedasticty: residuals are equal across the regression line
# Normality of residuals: residuals of the regression line are approximately normally distributed

# The dataset contains disease and physical indicators

# Task1: Explore the relationship between BMXBMI (BMI) and RIDAGEYR (age)
data <- read.csv('week7_data.csv')
summary(data)
data <- subset(data, !is.na(BMXBMI))
f <- ggplot(data)
# Check the distribution and linearity between 2 variables
f + geom_histogram(aes(RIDAGEYR))
f + geom_histogram(aes(BMXBMI))
f + geom_point(aes(RIDAGEYR, BMXBMI))

#Pre processing, remove outliers and confine the explore range to adults
data <- subset(data, data$BMXBMI < 100)
data_a <- subset(data, data$RIDAGEYR > 20)

#check data again
f <- ggplot(data_a)
f + geom_histogram(aes(RIDAGEYR))
f + geom_histogram(aes(BMXBMI))
f + geom_point(aes(RIDAGEYR, BMXBMI))

# Normal distribution test
qqnorm(data$RIDAGEYR)
qqline(data$RIDAGEYR)
qqnorm(data$BMXBMI)
qqline(data$BMXBMI)

# Spearman Correlation analysis ( do not requires normal distribution)
cor(data_a$RIDAGEYR, data_a$BMXBMI, method = "spearman")

# Linear regression
fit1 <- lm(data_a$BMXBMI~data_a$RIDAGEYR)
summary(fit1)
fit2 <- lm(data_a$BMXBMI~data_a$RIDAGEYR + I(data_a$RIDAGEYR ^ 2))
summary(fit2)
fit3 <- lm(data_a$BMXBMI~data_a$RIDAGEYR + I(data_a$RIDAGEYR^3))
summary(fit3)

# Check homoscedasticity and normality of residuals
glejser(fit3)
res1 = data.frame(residuals(fit1))
ggplot(res1) + geom_histogram(aes(residuals.fit1.))
qqnorm(res1$residuals.fit1.)
qqline(res1$residuals.fit1., col = 'red')


# Task2: Add other variables and see how they impact the linear model
# Add more independent variables
data_a_i = data_a[c('BMXBMI', 'RIDAGEYR', 'BMXWAIST', 'LBXGH', 'LBDTRSI', 'LBDLDLSI', 'gender_01', 'smoke_01')]
data_a_i = na.omit(data_a_i)

# Create histogram for each variable
# Create scatter plot between each variable and BMI to see linearity
ggplot(data_a_i) + geom_histogram(aes(smoke_01))
ggplot(data_a_i) + geom_point(aes(LBDTRSI, BMXBMI))

# Remove outliers
data_a_i <- subset(data_a_i, data_a_i$BMXWAIST<999)
data_a_i <- subset(data_a_i, data_a_i$LBXGH<13)

# Correlation coefficient analysis
# For normal distributed variables, perform Pearson correlation analysis
# For non-normal distributed variables, perform Spearman correlation analysis
cor(data_a_i$BMXWAIST, data_a_i$BMXBMI, method = "pearson")
cor(data_a_i$LBXGH, data_a_i$BMXBMI, method = "pearson")
cor(data_a_i$LBDTRSI, data_a_i$BMXBMI, method = "pearson")
cor(data_a_i$LBDLDLSI, data_a_i$BMXBMI, method = "pearson")
cor(data_a_i$gender_01, data_a_i$BMXBMI, method = "spearman")
cor(data_a_i$smoke_01, data_a_i$BMXBMI, method = "spearman")

# Fit the model
fit4 <- lm(data_a_i$BMXBMI ~ data_a_i$LBDTRSI + data_a_i$LBXGH + data_a_i$BMXWAIST + data_a_i$LBDLDLSI + data_a_i$RIDAGEYR + I(data_a_i$RIDAGEYR^3) + data_a_i$gender_01 + data_a_i$smoke_01)
summary(fit4)
# Fit the model without cubic term
fit5 <- lm(data_a_i$BMXBMI ~ data_a_i$LBDTRSI + data_a_i$LBXGH + data_a_i$BMXWAIST + data_a_i$LBDLDLSI + data_a_i$RIDAGEYR + data_a_i$gender_01 + data_a_i$smoke_01)
summary(fit5)
# Fit the model without cubic term and binary variables
fit6 <- lm(data_a_i$BMXBMI ~ data_a_i$LBDTRSI + data_a_i$LBXGH + data_a_i$BMXWAIST + data_a_i$LBDLDLSI + data_a_i$RIDAGEYR)
summary(fit6)

# Create coefficient matrix to exam multicollinearity
data_a_i_m = data_a[c('RIDAGEYR', 'BMXWAIST', 'LBXGH', 'LBDTRSI', 'LBDLDLSI', 'gender_01', 'smoke_01')]
data_a_i_m = na.omit(data_a_i_m)
data_a_i_m <- subset(data_a_i_m, data_a_i_m$BMXWAIST<999)
data_a_i_m <- subset(data_a_i_m, data_a_i_m$LBXGH<13)
cor(data_a_i_m, method = 'pearson')

# Check auto-correlation, homoscedasticity and normality of residuals

glejser(fit5)
dwtest(fit5)

res5 = data.frame(residuals(fit5))
ggplot(res5) + geom_histogram(aes(residuals.fit5.))
qqnorm(res5$residuals.fit5.)
qqline(res5$residuals.fit5., col = 'red')
shapiro.test(res5$residuals.fit5.)


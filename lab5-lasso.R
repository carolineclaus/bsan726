## Lab 5 -- LASSO
library(MASS)
library(dplyr)
library(glmnet)

## get data ready
data(Boston)

Boston$rad <- as.factor(Boston$rad)

Y <- Boston$medv
X <- dplyr::select(Boston, -medv)
X_std <- X %>% mutate(across(where(is.numeric), scale))

Xmat <- model.matrix(~., X_std)
head(Xmat)
Xmat <- Xmat[,-1] #We do not need the first column (intercept).

## Fitting lasso model
install.packages("glmnet")
library(glmnet)
lasso.fit<- glmnet(x=Xmat, y=Y)
# solution path
plot(lasso.fit, xvar = "lambda")
# coefficients
coef_lamb1 <- coef(lasso.fit, s=0.1)  # s specifies lambda value
coef_lamb2 <- coef(lasso.fit, s=0.5)
coef_lamb3 <- coef(lasso.fit, s=1)

# cross-validation
cv.lasso<- cv.glmnet(x=Xmat, y=Y)
plot(cv.lasso)
cv.lasso$lambda.min
cv.lasso$lambda.1se
# optimal model
coef.lamb1se <- coef(lasso.fit, s=cv.lasso$lambda.min)

## model assessment
## MSE
# first get predicted Y 
# note that we use the same predict() function, but the input is slightly different from the one we used for lm.
pred.lambmin <- as.vector(predict(lasso.fit, newx = Xmat, s=cv.lasso$lambda.min))
# then figure out degrees of freedom
df.lambmin <- nrow(Xmat)-sum(coef.lambmin!=0)
# then compute mse
rss.lambmin <- sum((pred.lambmin-Y)^2)
mse.lambmin <- rss.lambmin/df.lambmin

## Adjusted R^2
n <- nrow(X)
tss <- var(Y)*(n-1)
adjR2 <- 1-(n-1)/df.lambmin*rss.lambmin/tss

## CV score
ind.lambmin <- cv.lasso$index[1]
cv.lambmin <- cv.lasso$cvm[ind.lambmin]

### Exercise
# 1. For the other optimal model (lambda 1se), what are selected variables? Obtain MSE, adjusted R^2 and CV score.
# 2. Randomly split the data to training (80%) and testing (20%) and conduct the training-testing procedure. Obtain MSPE and oos-R^2.

### Exercise: explore group lasso using package "gglasso".

#---------------------------------------

## Ridge regression
ridge.fit<- glmnet(x=Xmat, y=Y, family = "gaussian", alpha = 0)
plot(ridge.fit, xvar = "lambda")

cv.ridge<- cv.glmnet(x=Xmat, y=Y, family = "gaussian", alpha = 0, nfolds = 10)
plot(cv.ridge)

pred.ridge.min<- predict(ridge.fit, newx = Xmat, s=cv.ridge$lambda.min)
rss.min <- sum((Y-pred.ridge.min)^2)
pred.ridge.1se<- predict(ridge.fit, newx = Xmat, s=cv.ridge$lambda.1se)
rss.1se <- sum((Y-pred.ridge.1se)^2)

#---------------------------------------

## High-dimensional regression
genedata<- read.csv("https://www.dropbox.com/s/idok8ibbcw5h5ql/gene_exp.csv?dl=1")
dim(genedata)
# OLS
lm.fit.high<- lm(Y~., data = genedata)
lm.fit.high
# LASSO
lasso.fit.high<- glmnet(x= as.matrix(genedata[,-1]), y=genedata[,1], family = "gaussian", alpha=1)
cv.lasso.fit.high<- cv.glmnet(x= as.matrix(genedata[,-1]), y=genedata[,1], family = "gaussian", alpha=1)
plot(cv.lasso.fit.high)

# obtain coefficient estimates
coef.high<- as.vector(coef(lasso.fit.high, s=cv.lasso.fit.high$lambda.1se))
coef.high[which(coef.high!=0)]

colnames(genedata[,-1])[which(coef.high!=0)]


## Bias-variance tradeoff (using ridge regression)





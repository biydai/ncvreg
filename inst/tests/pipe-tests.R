
# test PIPE function

#------------------------------------
# test linear case
#------------------------------------

# test 1
n <- 100
p <- 200
X <- matrix(rnorm(n*p),nrow = n)
y <- X %*% c(rep(0.5,10),rep(0,p - 10)) + rnorm(n)

data <- list(X = X, y = y)
cv.fit <- cv.ncvreg(X,y,penalty = "MCP",nfolds = 10,family = "gaussian")
fit <- ncvreg(data$X,data$y,penalty = "MCP",nfolds = 10,family = "gaussian")
pipe_result <- pipe(fit, cv.fit$lambda.min)
head(pipe_result$t,20)

# test 2
n <- 100
p <- 200
X <- matrix(rnorm(n*p),nrow = n)
y <- X %*% c(rep(0.5,10),rep(0,p - 10)) + rnorm(n)
data <- list(X = X, y = y)
cv.fit <- cv.ncvreg(X,y,penalty = "MCP",nfolds = 10,family = "gaussian")
fit <- ncvreg(data$X,data$y,penalty = "MCP",nfolds = 10,family = "gaussian", returnX = FALSE)
pipe_result <- pipe(fit, cv.fit$lambda.min)

# test 3
n <- 50
p <- 20
X <- matrix(rnorm(n*p), n, p)
b <- c(1, -1, 0.5, -0.5, rep(0.25, 3), rep(-0.25, 3), rep(0, 10))
y <- rnorm(n, X%*%b)

fit <- ncvreg(X, y, returnX=TRUE)
pipe(fit, lambda=0.2,sigmasq = 1)
fit <- ncvreg(X, y, returnX=TRUE, penalty='lasso')
pipe(fit, lambda=0.2,sigmasq = 1)

#------------------------------------
# test binomial case
#------------------------------------
set.seed(100)
n <- 100
p <- 20
X <- matrix(rnorm(n*p), n, p)
b <- c(1, -1, 0.5, -0.5, rep(0.25, 3), rep(-0.25, 3), rep(0, 10))
y <- rnorm(n, X%*%b) > 0
cv.fit <- cv.ncvreg(X,y,nfolds = 10,family = "binomial")
fit <- ncvreg(X, y, family='binomial', returnX=TRUE)
pipe(fit, lambda=cv.fit$lambda.min)$sigmasq

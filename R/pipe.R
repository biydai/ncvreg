#' Compute projection-based inference test statistics for all features
#' 
#' `pipe()` is a function that compute projection-based inference test statistics for all features
#' 
#' @param fit      A fitted `ncvreg` object.
#' @param lambda   The value of lambda at which inference should be carried out.
#' @param X,y      The design matrix and response used to fit the model; in most cases, it is not necessary to provide
#'   `X` and `y` as they are returned by `ncvreg`, but see the `returnX` argument in [ncvreg()].
#' @param sigmasq The estimated error variance for linear model. Only needed for gaussian case. If not supplied, the residual sums of squares will be computed
#'   
#' @return The object returns a data frame with one row per feature and six columns:
#' * `coef.model`: The coefficient estimate from the penalized regression fit
#' * `coef.pipe`: The pipe coefficient
#' * `SE`: Standard error of the pipe coefficient
#' * `t`: The pipe test statistic, computed based on coef.pipe and SE. It approximately follows a standard normal distribution under the null hypothesis that the
#'        feature is independent of the outcome.
#' * `p.value`: The pvalue computed based on the test statistic t
#' * `p.adjust`: The qvalue computed based on Benjamini Hochberg Procedure
#' 
#' @examples
#' # Linear regression
#' data(Prostate)
#' fit <- ncvreg(Prostate$X, Prostate$y)
#' pipe(fit, 0.1)
#' 
#' 
#' # Logistic regression
#' data(Heart)
#' X <- Heart$X
#' y <- Heart$y
#' fit <- ncvreg(X, y, family='binomial')
#' pipe(fit, 0.1)


pipe <- function(fit,# ncvreg object
                 lambda, # the model of interest, recommend: use cv
                 X = NULL,# design matrix
                 y = NULL,# outcome
                 sigmasq = NULL# this is an option for customized estimates
){ 
  # setting up
  stopifnot(class(fit)[1] == "ncvreg")
  
  # load data
  if (is.null(X) & is.null(fit$X)) {
    stop("This procedure requires X and y.  Either supply X and y, or fit the model using the option 'returnX = TRUE'")
  }
  
  tmp <- if (is.null(fit$X)) ncvreg(X, y, family=fit$family) else fit
  XX <- tmp$X
  yy <- tmp$y
  
  p <- ncol(XX)
  n <- nrow(XX)
  
  # load model
  beta <- coef(fit,lambda = lambda)
  
  S_hat <- which(beta[-1] != 0)
  N_hat <- which(beta[-1] == 0)
  
  # initialize
  beta_PIPE <- numeric(p)
  sigmasq_PIPE <- numeric(p)
  
  if(fit$family == "gaussian"){
    
    # compute sigmasq if not provided
    if(is.null(sigmasq)){
      sigmasq <- crossprod(yy - cbind(rep(1,n),XX) %*% beta)/(n - sum(beta[-1] != 0))
    }
    
    # compute pipe for features in the estimated support
    for (i in S_hat){
      S_hat_i <- beta[-1]!=0
      S_hat_i[i] <- FALSE
      #variance
      Xsi <- (XX[,S_hat_i,drop = FALSE])
      Qsi <- diag(n) - tcrossprod(qr.Q(qr(Xsi)))
      beta_PIPE[i] <- XX[,i] %*%(yy - Xsi %*% beta[-1][S_hat_i])/n 
      sigmasq_PIPE[i] <- sigmasq/(t(XX[,i])%*% Qsi %*% XX[,i])
    }
  
    # compute pipe for null features
    Xs <- XX[,S_hat]
    Qs <- diag(n) - tcrossprod(qr.Q(qr(Xs)))
  
    for(i in N_hat){
      beta_PIPE[i] <- XX[,i] %*%(yy - Xs %*% beta[-1][S_hat])/n 
      sigmasq_PIPE[i] <- sigmasq/(t(XX[,i])%*% Qs %*% XX[,i])
    }
  }  
  
  if(fit$family == "binomial"){
    
    # compute pseudo outcome
    pii <- predict(fit,XX,type = "response",lambda = lambda)
    A <- pii*(1-pii)
    W <- diag(A) # 
    v <- yy - pii # W and v change for different family of distribution
    Y_pse <- diag(1/A) %*% v + predict(fit,XX,type = "link",lambda = lambda)
    
    # compute weight matrix
    sqrtW <- sqrt(W)
    yy_pse <- sqrtW %*% Y_pse
    X_int <- cbind(rep(1,n),  XX)
    
    # compute pipe for features in the estimated support
    for (i in S_hat){
      S_hat_i <- beta[-1]!=0
      S_hat_i[i] <- FALSE
      
      Xs <- X_int[,c(TRUE,S_hat_i),drop = FALSE]
      XXs <- sqrtW %*% Xs
      xx <- crossprod(sqrtW,XX[,i])
      
      Qs <- diag(n) - tcrossprod(qr.Q(qr(XXs)))
      beta_PIPE[i] <- t(xx) %*% (yy - XXs %*% beta[c(TRUE,S_hat_i)])/crossprod(xx) 
      sigmasq_PIPE[i] <- 1/(t(xx) %*% Qs %*% xx)
    }
    
    # compute pipe for null features
    S_hat <- beta[-1]!=0
    Xs <- X_int[,c(TRUE,S_hat),drop = FALSE]
    XXs <- sqrtW %*% Xs
    Qs <- diag(n) - tcrossprod(qr.Q(qr(XXs)))
    
    for(i in N_hat){
      xx <- crossprod(sqrtW,XX[,i])
      beta_PIPE[i] <- t(xx) %*% (yy - XXs %*% beta[c(TRUE,S_hat)])/crossprod(xx) 
      sigmasq_PIPE[i] <- 1/(t(xx) %*% Qs %*% xx)
    }
    
  }
  
  t <- beta_PIPE/sqrt(sigmasq_PIPE)
  pvalue <-(1 - pnorm(abs(t)))*2
  qvalue <- p.adjust(pvalue,method = "BH")
  
  res <- data.frame(
      coef.model = beta[-1],
      coef.pipe = beta_PIPE,
      SE = sqrt(sigmasq_PIPE),
      t = t,
      p.value = pvalue,
      p.adjust = qvalue)
  
  res_reorder <- res[order(pvalue),]
  
  return(list(model = data.frame(lambda = lambda),
              results = res_reorder)
              
  )
}
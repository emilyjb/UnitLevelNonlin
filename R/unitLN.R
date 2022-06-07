#' Small area inference for the unit-level lognormal model
#'
#' @param yspos a numeric vector with the positive (not log transformed) response variables
#' @param Xs the matrix of covariates for sampled elements that does not contain an intercept
#' @param Xpop the matrix of covariates for the full population that does not contain an intercept
#' @param areafacpop a vector of area labels for the full population
#' @param areafacsamp a vector of area labels for the sample
#' @param sampindex the vector of sampled index values
#'
#' @return a list with predictions and MSE estimates
#' @export
#'
#' @examples
#' beta0 <- log(0.5)
#' beta1 <- 1
#' beta2 <- 2
#' D <- 60
#' Nis <- c(100, 200, 500)
#' Nis <- rep(Nis, each = D/length(Nis))
#' N <- sum(Nis)
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' areafacpop <- rep(1:D, Nis)
#' sigma2b <- 0.5
#' sigma2e <- 1
#' bi <- rnorm(D, mean = 0, sd = sqrt(sigma2b))
#' ei <- rnorm(N, mean = 0, sd = sqrt(sigma2e))
#' names(bi) <- as.character(1:D)
#' y <- exp(beta0 + beta1*x1 + beta2*x2 + bi[as.character(areafacpop)] + ei)
#' samplist <- sapply(1:D, function(i){sample( (1:N)[areafacpop == i],
#' size = 0.1*Nis[i], replace = FALSE)})
#' sampindex <- unlist(samplist)
#' areafacsamp <- areafacpop[sampindex]
#' ys <- y[sampindex]
#' Xpop <- cbind(x1, x2)
#' Xs <- Xpop[sampindex,]
#' popindex <- 1:N
#' unitLN(ys, Xs, Xpop,  areafacpop, areafacsamp, sampindex)
unitLN <- function(yspos, Xs, Xpop,  areafacpop, areafacsamp, sampindex){

  ys <- log(yspos)

  fitlog <- lme4::lmer(ys~Xs + (1|areafacsamp))
  sig2bhat <- lme4::VarCorr(fitlog)[[1]][1,1]
  sig2ehat <- summary(fitlog)$sigma^2
  betahat <- lme4::fixef(fitlog)

  nis <- tapply(areafacsamp, areafacsamp, length)
  Nis <- tapply(areafacpop, areafacpop, length)
  Gs <- stats::model.matrix(~ as.factor(areafacsamp) - 1)
  GN <- stats::model.matrix(~ as.factor(areafacpop)  -1)

  lbarsi <- tapply(ys, areafacsamp, mean)
  dbarsi <- t(Gs)%*%cbind(1, Xs)/as.vector(nis)

  gammais <- sig2bhat/(sig2bhat + sig2ehat/nis)
  btermis <- 0.5*(gammais/nis*sig2ehat + sig2ehat)
  sumnobs <- t(GN[-sampindex,])%*%(exp(Xpop[-sampindex,]%*%betahat[-1]))
  cis <- exp(betahat[1] + gammais*(lbarsi - as.vector(dbarsi%*%betahat))  + btermis)

  yhatimmse <- (as.vector(t(Gs)%*%exp(ys)) + cis*sumnobs[,1])/Nis

  ## Components for estimated MSE

  kappai <- exp(2*betahat[1] + 2*gammais^2*(sig2bhat + sig2ehat/nis) + gammais/nis*sig2ehat + sig2ehat)
  psi <- exp(gammais/nis*sig2ehat + sig2ehat) - exp(gammais/nis*sig2ehat)
  xi <- exp(gammais/nis*sig2ehat) - 1
  sumnobs2 <- t(GN[-sampindex,])%*%(exp(2*Xpop[-sampindex,]%*%betahat[-1]))

  yhatijr <- (GN[-sampindex,]%*%cis)*(exp(Xpop[-sampindex,]%*%betahat[-1]))
  bij1 <- 1-GN[-sampindex,]%*%gammais
  bij2 <- Xpop[-sampindex, ] - GN[-sampindex,]%*%(dbarsi[,-1]*as.vector(gammais))
  ais <- yhatijr[,1]*cbind(bij1, bij2)
  air <- t(GN[-sampindex,])%*%ais

  Is <- 0.5*apply(apply(GN[sampindex,], 2, vcovraneff2, sig2bhat, sig2ehat ),1,sum)
  Imat <- matrix(Is,2,2)
  vcovsig <- solve(Imat)

  fi1 <- 1/(sig2bhat + sig2ehat/nis)*(1-gammais)*(lbarsi - as.vector(dbarsi%*%betahat)) + 0.5*(1-gammais)^2
  fi2 <- -gammais/(sig2bhat + sig2ehat/nis)/nis*(lbarsi - as.vector(dbarsi%*%betahat)) + 0.5*gammais^2/nis + 0.5
  yhatir <- t(GN[-sampindex,])%*%yhatijr
  mse.sig <- diag(cbind(fi1, fi2)%*%vcovsig%*%t(cbind(fi1, fi2)))/Nis^2*as.vector(yhatir)^2
  M <- as.matrix(stats::vcov(fitlog))
  mse.beta <- apply(air, 1, function(v){ t(v)%*%M%*%v })/Nis^2

  MSEhatln <- kappai/Nis^2*(as.vector(sumnobs)^2*xi + as.vector(sumnobs2)*psi)

  list( fixpar = list(betahat, sig2bhat, sig2ehat), pred = yhatimmse , msehat = MSEhatln)
}


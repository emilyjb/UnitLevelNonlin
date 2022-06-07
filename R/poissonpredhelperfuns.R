
lpgd <- function(par, ysd, xsd){
  alpha <- par[1]; beta <- par[2]; gamma <- par[3]
  lambdaij <- exp(xsd*gamma)
  lambdaiprod <- prod(lambdaij^ysd)
  lambdaisum <- sum(lambdaij)
  alpha*log(beta) - lgamma(alpha) + sum(ysd*log(lambdaij)) + lgamma(sum(ysd) + alpha) - (sum(ysd) + alpha)*log(( beta + lambdaisum))
}

lpg <- function(par, ys, xs, areafacsamp){
  areas <- unique(areafacsamp)
  ls <- sapply(areas, function(d){   lpgd(par, ys[areafacsamp == d], xs[areafacsamp == d])})
  -sum((ls))
}



InitvalsPoisson <- function(ys, xs, areafacsamp  ){
  D <- length(table(areafacsamp))
  poissonglm <- stats:glm(ys~xs + as.factor(areafacsamp) - 1, family = poisson(link = "log"))
  lambdahatinit <- exp(xs%*%poissonglm$coef[1:ncol(xs)])
  lambdaidotinit <- tapply(lambdahatinit, areafacsamp, sum)
  thetahatinit <- mean(tapply(ys, areafacsamp, sum))/mean(lambdaidotinit)
  yidot <- tapply(ys, areafacsamp, sum)
  alphabybeta2init <- (mean( (yidot - lambdaidotinit*thetahatinit)^2) - mean(lambdaidotinit*thetahatinit))/mean(lambdaidotinit^2)
  betahatinit <- thetahatinit/alphabybeta2init
  alphahatinit <- thetahatinit*betahatinit

  checkdelta <- betahatinit > 0 & alphahatinit > 0

  if(!checkdelta ){
    uhatinit <- exp(summary(poissonglm)$coef[,"Estimate"][-c(1:ncol(xs))])
    muhatinit <- mean(uhatinit)
    alphahatinit <- stats:var(uhatinit)/muhatinit^2
    betahatinit <- alphahatinit/muhatinit
  }

  gammahatinit <- summary(poissonglm)$coef[,"Estimate"][1:ncol(xs)]

  c(alphahatinit, betahatinit, gammahatinit)

}


getEstPoisson <- function(parinit, ys, xs, areafacsamp){
  D <- length(table(areafacsamp))
  parestopt <- list(par = NULL, convergence = 1, hessian = NULL)
  try(parestopt <- stats:optim( lpg , par = parinit,  ys = ys, xs = xs, areafacsamp = areafacsamp, hessian = TRUE))
  parest <- parestopt$par
  H <- parestopt$hessian

  if(parestopt$convergence != 0){

    parest <- parinit

    #parest <- c(alphahatinit, betahatinit, gammahatinit)
  }

  list(parest = parest, hessian = H,  converges = parestopt$convergence )

}


PredG1G2Poisson <- function(parest, ys, xs, areafacsamp, xN, areafacpop, sampindex){
  D <- length(table(areafacsamp))
  alphahat <- parest[1]; betahat <- parest[2]; gammahat <- parest[-c(1,2)]

  yisum <- tapply(ys, areafacsamp, sum)
  lambdaijsamp <- exp(xs%*%gammahat)
  lambdaisum <- tapply(lambdaijsamp, areafacsamp , sum)

  uhat <- (yisum + parest[1])/(betahat + lambdaisum)
  names(uhat) <- 1:D

  yijhatpop <- exp(xN%*%gammahat)*uhat[areafacpop]
  yijhatpop[sampindex] <- ys

  ybarNihatEBfixed <- tapply(yijhatpop, areafacpop, mean)

  Nis <- tapply(areafacpop, areafacpop, length)
  ##### Compute g1fixed:
  d <- 0
  g1hatfixed  <- c()

  repeat{

    d <- d + 1
    xr <- xN[-sampindex,]
    areafacpopr <- areafacpop[-sampindex]

    xrd <- xr[areafacpopr == d,]

    fac1d <- (yisum[d] + alphahat)/(betahat + lambdaisum[d])
    term1d <- 1/Nis[d]^2*sum(exp(xrd%*%gammahat))*fac1d
    fac2d <- (yisum[d] + alphahat)/(betahat + lambdaisum[d])^2

    matx1 <- kronecker(t(rep(1, nrow(xrd))), as.vector(xrd%*%gammahat))
    matx2 <- t(matx1)

    matx <- matx1 + matx2

    term2d <- fac2d/Nis[d]^2*sum(exp(matx))

    g1hatfixed  <- c( g1hatfixed  , term1d + term2d)

    if(d == D){break}

  }

  g2mean <- NULL
  if(!is.null(H)){

    ### Taylor linearization variance estimation:
    vhatdeltahat <- solve(H)
    d <- 0
    g2mean <- c()
    repeat{
      d <- d + 1
      xNd <- xN[-sampindex[areafacpop[-sampindex] == d],]
      der1d1 <- apply(as.vector(exp(xNd%*%gammahat))*xNd, 2, sum)*(yisum[d] + alphahat)/(betahat + lambdaisum[d])/Nis[d]
      der1d2 <- -(yisum[d] + alphahat)/(betahat + lambdaisum[d])^2*sum(exp(xNd%*%gammahat))*apply(exp(xs[areafacsamp==d,]%*%gammahat)*xs[areafacsamp==d,], 2, sum)/Nis[d]
      der1d <- der1d1 + der1d2
      der2 <- -(yisum[d] + alphahat)/(betahat + lambdaisum[d])^2*sum(exp(xNd%*%gammahat))/Nis[d]
      der3 <- sum(exp(xNd%*%gammahat))/(betahat + lambdaisum[d])/Nis[d]
      g2d <- (t(c(der3, der2,der1d))%*%vhatdeltahat%*%c(der3, der2, der1d))[1,1]
      g2mean <- c(g2mean, g2d)
      if(d == D){break}
    }

  }

  list(ybarNihatEBfixed, g1hatfixed, g2mean)

}


#' Small area inference for a unit-level gamma-Poisson model
#'
#' @param ys a numeric vector containing the collected counts
#' @param xs the matrix of covariates for sampled elements that does not contain an intercept
#' @param areafacsamp a vector of area labels for the sample
#' @param xN the matrix of covariates for population elements that does not contain an intercept
#' @param areafacpop a vector of area labels for the population
#' @param sampindex indexes of sampled elements
#'
#' @return a list with four components: parest = estimators of fixed parameters, pred = predictors, g1hat = estimate of leading term in MSE, g2hat = estimate of second term in MSE
#' @export
#'
#' @examples
#' #### Simulate population:
#' N <- 10000
#' mux <- 0.5; sdx <- 1
#' D <- 100
#' Nis <- rep(100,D)
#' f <- 0.05
#' nis <- Nis*f
#' xN1 <- rnorm(N, mean = mux, sd = sdx)
#' xN2 <- rnorm(N, mean = mux, sd = sdx)
#' beta1 <- 2
#' muxN <- exp(beta1*xN1/4 + beta1*xN2/4)
#' uD <- rgamma(D, shape = 5, rate = 2)
#' names(uD) <- 1:D
#' areafacpop <- rep(1:D, Nis)
#' meanyN <- muxN*uD[as.character(areafacpop)]
#' yN <- rpois(N, lambda = meanyN)
#' indexpop <- 1:N
#' ybarNipop <- tapply(yN, areafacpop, mean)
#' sampindex <-  as.vector(sapply(1:D, function(i){ sample(indexpop[areafacpop == i],
#'      size = nis[i], replace = FALSE) }) )
#' xN <- cbind(xN1, xN2)
#' ys <- yN[sampindex]
#' xs <- xN[sampindex,]
#' areafacsamp <- areafacpop[sampindex]
SAEPoissonMean <- function( ys, xs, areafacsamp, xN, areafacpop, sampindex){
  parinit <- InitvalsPoisson(ys, xs, areafacsamp  )
  parestPoisson <- getEstPoisson(parinit, ys, xs, areafacsamp)
  predg1g2hat <- PredG1G2Poisson(parestPoisson, ys, xs, areafacsamp, xN, areafacpop, sampindex)
  list(parest = parestPoisson, pred = predg1g2hat[[1]], g1hat = predg1g2hat[[2]], g2hat = predg1g2hat[[3]])
}


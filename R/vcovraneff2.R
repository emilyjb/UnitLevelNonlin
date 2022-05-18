vcovraneff2 <- function(x, sig2u, sig2e){
  ni <- sum(x)
  alphai <- sig2e + ni*sig2u
  Ivv <- ni^2/alphai^2
  Iee <- (ni-1)/sig2e^2 + 1/alphai^2
  Ive <- ni/alphai^2
  c(Ivv, Ive, Ive, Iee)
}

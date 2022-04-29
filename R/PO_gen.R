# Generate data from PO

# lp: beta * Z
# inv.m: inverse of m(t)
PO.gen = function(lp, inv.m)
{
  return(inv.m(logit(runif(length(lp))) - lp))
}


#' generate data for the example
#' @param n Sample size
#' @param beta True parameter, p-dimensional vector
#' @param Z.gen Function to generate n by p covariate matrix
#' @param inv.m Inverse of baseline hazard function
#' @param C.gen Function to generate n-dimensional censoring time
#' @export
PO.sim = function(n, beta,
                  Z.gen = function(n, p) matrix(rnorm(n*p), n, p),
                  inv.m = exp,
                  C.gen = function(n) runif(n, 0, 5))
{
  p = length(beta)

  Zn = Z.gen(n,p)
  lp = drop(Zn %*% beta)

  Tn = PO.gen(lp, inv.m)
  Cn = C.gen(n)

  dn = as.numeric(Tn <= Cn)
  Xn = pmin(Tn, Cn)

  return(list(delta = dn, X = Xn,
              Z = Zn, C = Cn))
}

# Utility functions
logit = function(x)
{
  log(1/(1/x-1))
}

expit = function(x)
{
  1/(1+exp(-x))
}

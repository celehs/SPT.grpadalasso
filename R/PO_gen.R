# Generate data from PO

# lp: beta * Z
# inv.m: inverse of m(t)
PO.gen = function(lp, inv.m)
{
  return(inv.m(logit(runif(length(lp))) - lp))
}


#' generate data for the example
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

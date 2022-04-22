#
# xi.CWY = function(x)
# {
#   pos = x>0
#   neg = x<0
#
#   out = rep(0.5, length(x))
#   eneg = exp(x[neg])
#   out[neg] = (eneg*x[neg]+1-eneg)/(1-eneg)^2
#   inv.epos = exp(-x[pos])
#   out[pos] = inv.epos*(x[pos]+inv.epos-1)/(inv.epos-1)^2
#
#   return(out)
# }

xi.CWY = function(x)
{
  ex = exp(-abs(x))
  out = (ex*(x-1) + ex^(2*(x>0)))/(1-ex)^2
  out[is.nan(out)] = 0.5
  return(out)
}

dxi.CWY = function(x)
{
  x = -abs(x)
  ex = exp(x)
  out = (x*ex*(ex+1) + 2*ex*(1-ex))/(1-ex)^3
  out[is.nan(out)] = -1/6
  return(out)
}

# x=seq(-5,5,0.1)
# plot(x,xi.deriv.fun(x), type = 'l', col =2 )
# lines(x,dxi.CWY(x), lty=2, col =3 )

#
# surv = Surv(dat$X, dat$delta)
# Z = dat$Z
# maxit = 1000
# tol = 1e-6
# GC = 1-ecdf(dat$C)(dat$X)
# beta = rep(0, ncol(Z))

PO.ipcw = function(surv, Z, GC,
                   beta, maxit = 1000, tol = 1e-6)
{
  n = length(surv)
  p = ncol(Z)

  # Initial value
  if(missing(beta))
  {
    beta = rep(0, p)
  }
  # fast calculation
  X.order = order(surv[,1], decreasing = TRUE)

  Ij = rep(0,n)
  Zj = matrix(0,n,p)
  j = 1
  for(i in 1:n)
  {
    while(surv[X.order[i],1] < surv[X.order[j],1])
    {
      Ij[X.order[j]] = Ij[X.order[n]]
      Zj[X.order[j],] = Zj[X.order[n],]
      j = j+1
    }

    Ij[X.order[n]] = Ij[X.order[n]] + 1
    Zj[X.order[n],] = Zj[X.order[n],] + Z[X.order[i],]
  }

  resp.fast = apply((Zj - Z*Ij)*surv[,2]/GC^2, 2, sum)/n^2

  i.tri = rep(2:n, 2:n-1)
  j.tri = unlist(sapply(2:n-1,function(x) 1:x))
  up.pos = (i.tri-1)*n + j.tri
  low.pos = (j.tri-1)*n + i.tri
  rep1n = rep(1,n)
  xi.mat = matrix(0.5,n,n)
  dxi.mat = matrix(-1/6,n,n)


  for(iter in 1:maxit)
  {
    # print(beta)
    lp = drop(Z %*% beta)
    xi.mat[up.pos] = xi.CWY(lp[i.tri] - lp[j.tri])
    xi.mat[low.pos] = 1-xi.mat[up.pos]

    xi.Sj = drop(xi.mat %*% rep1n)
    pred.fast = drop((n-2*xi.Sj) %*% Z)/n^2

    U = resp.fast - pred.fast
    print('u:')
    print(U)
    if(max(abs(U)) < tol)
    {
      break
    }

    dxi.mat[up.pos] = dxi.CWY(lp[i.tri] - lp[j.tri])
    dxi.mat[low.pos] = dxi.mat[up.pos]

    dxi.Sj = drop(dxi.mat %*% rep1n)
    dxi.Zj = dxi.mat %*% Z

    Hess = 2*t(Z)%*%(dxi.Sj*Z-dxi.Zj)/n^2

    beta = beta + solve(Hess,resp.fast - pred.fast)
  }

  if(iter > maxit)
  {
    stop("Fail to converge at max iteration.")
  }

  return(beta)
}
#
# tseq = dat$X
# surv = Surv(dat$X, dat$delta)
# Z = dat$Z
# Gt = GX
# beta = beta.CWY
# maxit = 1000
# tol = 1e-6
# max.move = 1

PO.base.ipcw = function(tseq, surv, Z,  Gt,
                    beta, maxit = 1000, tol = 1e-8, max.move = 1)
{
  k = length(tseq)
  n = length(surv)
  t.order = order(tseq, decreasing = T)
  X.order = order(surv[,1], decreasing = T)
  elp = exp(drop(Z %*% beta))

  X.pos = 1
  cumI = 0
  a = 0
  ht = rep(0,k)
  for (i in 1:k)
  {
    while(X.pos <= n)
    {
      if(surv[X.order[X.pos],1] <= tseq[t.order[i]])
      {
        break
      }
      cumI = cumI + 1
      X.pos = X.pos + 1
    }
    resp = cumI/n/Gt[t.order[i]]
    if(resp >= 1)
    {
      ht[t.order[i]] = -Inf
      next
    }
    if(resp <= 0)
    {
      ht[t.order[i]] = Inf
      next
    }

    for(iter in 1:maxit)
    {
      eq = resp - mean(1/(1+elp*a))
      if(abs(eq)<tol)
        break

      deq = mean(elp/(1+elp*a)^2)
      update = - eq/deq
      update = update* min(1, 1/abs(update))
      a = a + update
      while(a < 0)
      {
        update = update/2
        a = a - update
      }
    }

    if(iter > maxit)
    {
      stop("Fail to converge at max iteration.")
    }

    ht[t.order[i]] = log(a)

  }
  return(ht)
}


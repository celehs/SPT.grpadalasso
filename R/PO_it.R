require(survival)
#
# surv = Surv(dat$X, dat$delta)
# Z = dat$Z
# df = 10
# maxit = 1
# tol = 1e-6
# max.move = 1
# glm.maxit = 1
# lam.type = "breslow"
# beta = beta.CWY


PO.it = function(surv, Z,  df = round(sqrt(nrow(Z))),
                 bs.knots,
                 beta, lam, GX, a,
                 maxit = 1000, tol = 1e-6,
                 glm.maxit = 1, max.move = 1,
                 order = 1,
                 lam.type = c("newton","coord-joint", "coord-loop","breslow"))
{
  n = length(surv)
  p = ncol(Z)

  if(missing(bs.knots))
  {
    bs.knots = quantile(surv[surv[,2]==1,1], probs = seq(0,1, length.out = df+2))
  }else{
    df = length(bs.knots)-2
  }
  interior = (surv[,1] >= bs.knots[2]) & (surv[,1] <= bs.knots[df+1])

  # Calculate the bases
  if(order == 1){
    bs.x = bs.1(surv[,1],bs.knots, df)
  }

  # Initial value
  if(missing(beta))
  {
    beta = rep(0, p)
  }
  if(missing(lam))
  {
    if(missing(GX))
    {
      Gfit = survfit(Surv(surv[,1], 1-surv[,2])~1)
      GX = stepfun(Gfit$time, c(1,Gfit$surv))(surv[,1])
    }
    mX = PO.base.ipcw(surv[,1], surv, Z, GX, beta)
    lam.init = coef(lm(mX ~ bs.x$Ibs,
                       subset = interior))
    a = lam.init[1]
    lam = pmax(0,lam.init[-1])
  }
  mX = drop(bs.x$Ibs %*% lam)
  lp =  drop(Z %*% beta)

  # Pseudo data
  d.pseudo = c(rep(0:1, c(n,sum(surv[,2]))))
  pseudo.pos = c(1:n , which(surv[,2]==1))
  Z.pseudo = Z[pseudo.pos,]

  if(missing(a))
  {
    a = coef(glm(d.pseudo ~ 1, family = binomial,
                 offset = (mX+lp)[pseudo.pos]))
  }
  lp = lp + a

  # Iterative algorithm
  lam.num = drop(surv[,2] %*% bs.x$bs)
  iter = 0

  a.old = a
  beta.old = beta
  best.loglik = -Inf
  cont.beta = TRUE
  while(cont.beta)
  {
    # print(beta[1:5])
    cont.lam = TRUE
    while (cont.lam)
    {
      # print(lam[1:5])

      if(lam.type[1] == "breslow")
      {
        piX = expit(lp+mX)
        lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs)
        update = lam.num/lam.denom - lam
        lam.new = lam + update* pmin(max.move/sqrt(sum(update^2)), 1)
      }else if(lam.type[1] == "coord-joint")
      {
        piX = expit(lp+mX)
        lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs)
        lam.hess = drop(((1+surv[,2])*piX*(1-piX))%*% bs.x$Ibs^2)
        update = -(lam - lam.num/lam.denom)/(1+lam.num*lam.hess/lam.denom^2)
        lam.new = lam + update* pmin(max.move/sqrt(sum(update^2)), 1)
      }else if(lam.type[1] == "coord-loop")
      {
        lam.new = lam
        for (j in 1:df)
        {
          piX = expit(lp+mX)
          lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs[,j])
          lam.hess = drop(((1+surv[,2])*piX*(1-piX))%*% bs.x$Ibs[,j]^2)
          update = -(lam[j] - lam.num[j]/lam.denom)/(1+lam.num[j]*lam.hess/lam.denom^2)
          lam.new[j] = lam[j] + update* pmin(max.move/abs(update), 1)
          mX = mX + bs.x$Ibs[,j]*(lam.new[j]-lam[j])
        }
      }else if(lam.type[1] == "newton")
      {
        piX = expit(lp+mX)
        lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs)
        lam.hess = (lam.num/lam.denom^2)*(t(bs.x$Ibs) %*% (((1+surv[,2])*piX*(1-piX)) * bs.x$Ibs))
        diag(lam.hess) = 1+diag(lam.hess)
        update = - solve(lam.hess , lam - lam.num/lam.denom)
        lam.new = lam + update* pmin(max.move/sqrt(sum(update^2)), 1)
        # print(lam.num[1:5])
        # print(lam.denom[1:5])
        # print(lam.hess[1:5,1:5])
        # print(update[1:5])
      }

      iter = iter + 1
      if(iter > maxit)
      {
        stop("Fail to converge at max iteration.")
      }
      cont.lam = max(abs(lam.new-lam)) > tol
      lam = lam.new
      mX = drop(bs.x$Ibs %*% lam)
    }

    new.loglik = PO.bs.loglik(surv[,2], Z, bs.x,
                              c(a,beta,lam))
    # print(c(new.loglik,best.loglik, new.loglik < best.loglik))
    if(new.loglik < best.loglik)
    {
      a = a.old + (a-a.old)/2
      beta = beta.old + (beta-beta.old)/2
      lp = a + drop(Z %*% beta)

      if(max(abs(beta-beta.old))<tol)
        break
      next
    }
    best.loglik = new.loglik
    # X.order = order(surv[,1])
    # plot(surv[X.order,1], mX[X.order], type = 'l')

    a.tmp = coef(glm(d.pseudo ~ 1, family = binomial,
                     offset = (mX+lp-a)[pseudo.pos]))
    defaultW <- getOption("warn")
    options(warn = -1)
    tmp.fit = glm(d.pseudo ~ Z.pseudo, family = binomial,
                 offset = mX[pseudo.pos]
                 ,control = glm.control(maxit = glm.maxit)
                 ,start = c(a.tmp, beta)
                 )
    options(warn = defaultW)
    update = coef(tmp.fit)-c(a.tmp,beta)
    rescale =  min(1, max.move/sqrt(sum(update^2)))
    iter = iter + tmp.fit$iter

    beta = beta + update[-1]*rescale
    a = a.tmp+update[1]
    lp = a + drop(Z %*% beta)
    cont.beta = max(abs(update)) > tol
    a.old = a
    beta.old = beta

    # print(best.loglik)

    if(iter > maxit)
    {
      stop("Fail to converge at max iteration.")
    }
  }
  if(max(abs(update[-1]))>tol*100)
  {
    warning("Algorithm may not converge.")
  }

  return(list(beta = beta, a = a, lam = lam,
              bs.x = bs.x, iter = iter))

}
#
# delta = Delta
# bs.x = fit.bs$bs.x
# abetalam = fit.adaglasso$coefficients

PO.bs.loglik = function(delta, Z, bs.x,
                        abetalam)
{
  if(is.null(dim(abetalam)))
  {
    abetalam = matrix(abetalam)
  }
  n = length(delta)
  lpX = cbind(1,Z,bs.x$Ibs) %*% abetalam
  dmX = bs.x$bs %*% abetalam[1+ncol(Z)+1:ncol(bs.x$bs),,drop = F]
  dmX[delta==0,] = 1
  dmX[apply(bs.x$bs==0, 1, all),] = 1

  out = apply(delta*(lpX + log(dmX)) - (1+delta)*log(1+exp(lpX)),2,mean)
  nan.pos = which(is.nan(out))
  if(any(nan.pos))
    out[nan.pos] = -Inf

  out
}

# surv = Surv(SX, Delta)
# bs.x = fit.bs$bs.x
# beta = fit.adaglasso$coefficients[1+1:ncol(Z),]
# lam = fit.bs$lam
# maxit = 1000
# tol = 1e-6
# glm.maxit = 1
# max.move = 1
# lam.type = "newton"

PO.bs.ploglik = function(surv, Z, bs.x,
                         beta, lam, GX,
                         maxit = 1000, tol = 1e-6,
                         glm.maxit = 1, max.move = 1,
                         lam.type = c("newton","coord-joint", "coord-loop","breslow"))
{
  if(is.null(dim(beta)))
  {
    beta = matrix(beta)
  }

  n = length(surv)
  tmpbeta = beta[,ncol(beta)]
  lp = drop(Z %*% tmpbeta)

  d.pseudo = c(rep(0:1, c(n,sum(surv[,2]))))
  pseudo.pos = c(1:n , which(surv[,2]==1))
  Z.pseudo = Z[pseudo.pos,]

  if(missing(lam))
  {
    if(missing(GX))
    {
      Gfit = survfit(Surv(surv[,1], 1-surv[,2])~1)
      GX = stepfun(Gfit$time, c(1,Gfit$surv))(surv[,1])
    }
    mX = PO.base.ipcw(surv[,1], surv, Z, GX, tmpbeta)
    lam.init = coef(lm(mX ~ bs.x$Ibs,
                       subset = interior))
    a = lam.init[1]
    lam = lam.init[-1]
  }else{
    mX = drop(bs.x$Ibs %*% lam)
    a = coef(glm(d.pseudo ~ 1, family = binomial,
                 offset = (mX+lp)[pseudo.pos]))
  }

  out.a = out.loglik = rep(NA, ncol(beta))
  out.lam = matrix(0, length(lam), ncol(beta))
  for (i in ncol(beta):1)
  {
    lp = a + drop(Z %*% beta[,i])

    # Iterative algorithm
    lam.num = drop(surv[,2] %*% bs.x$bs)
    iter = 0

    cont.lam = TRUE
    while (cont.lam)
    {
      # print(lam[1:5])

      if(lam.type[1] == "breslow")
      {
        piX = expit(lp+mX)
        lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs)
        update = lam.num/lam.denom - lam
        lam.new = lam + update* pmin(max.move/sqrt(sum(update^2)), 1)
      }else if(lam.type[1] == "coord-joint")
      {
        piX = expit(lp+mX)
        lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs)
        lam.hess = drop(((1+surv[,2])*piX*(1-piX))%*% bs.x$Ibs^2)
        update = -(lam - lam.num/lam.denom)/(1+lam.num*lam.hess/lam.denom^2)
        lam.new = lam + update* pmin(max.move/sqrt(sum(update^2)), 1)
      }else if(lam.type[1] == "coord-loop")
      {
        lam.new = lam
        for (j in 1:df)
        {
          piX = expit(lp+mX)
          lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs[,j])
          lam.hess = drop(((1+surv[,2])*piX*(1-piX))%*% bs.x$Ibs[,j]^2)
          update = -(lam[j] - lam.num[j]/lam.denom)/(1+lam.num[j]*lam.hess/lam.denom^2)
          lam.new[j] = lam[j] + update* pmin(max.move/abs(update), 1)
          mX = mX + bs.x$Ibs[,j]*(lam.new[j]-lam[j])
        }
      }else if(lam.type[1] == "newton")
      {
        piX = expit(lp+mX)
        lam.denom = drop(((1+surv[,2])*piX-surv[,2])%*% bs.x$Ibs)
        lam.hess = (lam.num/lam.denom^2)*(t(bs.x$Ibs) %*% (((1+surv[,2])*piX*(1-piX)) * bs.x$Ibs))
        diag(lam.hess) = 1+diag(lam.hess)
        update = - solve(lam.hess , lam - lam.num/lam.denom)
        lam.new = lam + update* pmin(max.move/sqrt(sum(update^2)), 1)
        # print(lam.num[1:5])
        # print(lam.denom[1:5])
        # print(lam.hess[1:5,1:5])
        # print(update[1:5])
      }

      iter = iter + 1
      if(iter > maxit)
      {
        stop("Fail to converge at max iteration.")
      }
      cont.lam = max(abs(lam.new-lam)) > tol
      lam = lam.new
      mX = drop(bs.x$Ibs %*% lam)
      a.new = coef(glm(d.pseudo ~ 1, family = binomial,
                   offset = (mX+lp-a)[pseudo.pos]))
      lp = lp - a + a.new
      a = a.new
    }

    dmX = bs.x$bs %*% lam
    dmX[surv[,2]==0] = 1
    dmX[apply(bs.x$bs==0, 1, all)] = 1

    out.a[i] = a
    out.lam[,i] = lam
    out.loglik[i] = mean(surv[,2]*(lp+mX + log(dmX)) -
                    (1+surv[,2])*log(1+exp(lp+mX)))
  }

  return(list(a = out.a, lam = out.lam, loglik = out.loglik))
}

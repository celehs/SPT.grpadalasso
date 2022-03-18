PO.glasso.fit = function (delta, X, Z, C, df  , order = 1,
               method = c('glasso','glasso-PLH'),control)
{
  method <- match.arg(method)
  if (missing(C)){
    KMfit = summary(survfit(Surv(X,delta)~1), time = sort(unique(X)))
    GX = stepfun(KMfit$time, c(1,KMfit$surv))(X)
  }
  else{
    GX = 1-ecdf(C)(X-1e-8)
  }
  nn = length(delta)

if(method == 'glasso')
{
  if(missing(df)){
    df = as.integer(nn**(1/3))
  }
  grp = seq(1,ncol(Z))
  #GX = 1-ecdf(C)(X-1e-8)
  beta.CWY = PO.ipcw(Surv(X, delta), Z,GX, maxit = control$ipcw.maxit,tol = control$ipcw.tol)
  fit.bs = PO.it(Surv(X, delta), Z, df = df,
                 beta = beta.CWY,
                 # lam = rep(0,10),
                 GX = GX,
                 glm.maxit = control$it.maxit, lam.type = control$it.lam.type)


  fit.glasso = PO.glasso(delta, Z, fit.bs$bs.x,
                         grp,
                         fit.bs$beta, fit.bs$lam, fit.bs$a,
                         ada = control$ada,
                         lambda.min = control$lambda.min, nlambda = control$nlambda)


  glasso.loglik = PO.bs.loglik(delta, Z, fit.bs$bs.x,
                               fit.glasso$coefficients)*nn
  glasso.aic = 2*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.loglik
  glasso.bic = log(nn)*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.loglik

  beta.aic = fit.glasso$coefficients[1+0:length(beta.CWY),
                                     which.min(glasso.aic)]
  beta.bic = fit.glasso$coefficients[1+0:length(beta.CWY),
                                     which.min(glasso.bic)]
  beta = list(beta.aic,beta.bic)
  names(beta) = c("beta.aic","beta.bic")
}
else if(method == 'glasso-PLH')
{
  if(missing(df)){
    df = as.integer(nn**(1/3))
  }
  grp = seq(1,ncol(Z))
  #GX = 1-ecdf(C)(X-1e-8)
  grp = c(1,2,1, 2+rep(1:5,each = 2))
  beta.CWY = PO.ipcw(Surv(X, delta), Z,GX, maxit = control$ipcw.maxit,tol = control$ipcw.tol)
  fit.bs = PO.it(Surv(X, delta), Z, df = df,
                 beta = beta.CWY,
                 # lam = rep(0,10),
                 GX = GX,
                 glm.maxit = control$it.maxit, lam.type = control$it.lam.type)
  fit.glasso = PO.glasso(delta, Z, fit.bs$bs.x,
                         grp,
                         fit.bs$beta, fit.bs$lam, fit.bs$a)
  glasso.refit = PO.bs.ploglik(Surv(X, delta), Z, fit.bs$bs.x,
                               fit.glasso$coefficients[1+1:ncol(Z),],
                               fit.bs$lam, GX,
                               maxit = control$bs.ploglik.maxit, tol = control$bs.ploglik.tol,
                               glm.maxit = control$bs.ploglik.glm.maxit, max.move = control$bs.ploglik.max.move,
                               lam.type = control$bs.ploglik.lam.type
  )
  glasso.ploglik = glasso.refit$loglik*nn


  glasso.paic = 2*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.ploglik
  glasso.pbic = log(nn)*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.ploglik

  beta.paic = fit.glasso$coefficients[1+0:length(beta.CWY),
                                      which.min(glasso.paic)]
  beta.pbic = fit.glasso$coefficients[1+0:length(beta.CWY),
                                      which.min(glasso.pbic)]
  beta = list(beta.paic,beta.pbic)
  names(beta) = c("beta.aic","beta.bic")
}
  res = list()
  res$coefficients = beta
  return(res)
}

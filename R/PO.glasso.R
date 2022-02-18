PO.glasso = function (delta, X, Z, C, df  , order = 1,
               method = c('U-method','B-spline','NPMLE'),control)
{
  method <- match.arg(method)
  if (missing(C)){
    print('generate GX')
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
  #GX = 1-ecdf(C)(X-1e-8)
  beta.CWY = PO.ipcw(Surv(X, delta), Z,GX)
  fit.bs = PO.it(Surv(X, delta), Z, df = df,
                 beta = beta.CWY,
                 # lam = rep(0,10),
                 GX = GX,
                 glm.maxit = 1, lam.type = "newton")
  grp = c(1,2,1, 2+rep(1:5,each = 2))

  fit.glasso = PO.glasso(delta, Z, fit.bs$bs.x,
                         grp,
                         fit.bs$beta, fit.bs$lam, fit.bs$a)


  glasso.loglik = PO.bs.loglik(delta, Z, fit.bs$bs.x,
                               fit.glasso$coefficients)*nn
  glasso.aic = 2*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.loglik
  glasso.bic = log(nn)*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.loglik

  beta.aic = fit.glasso$coefficients[1+0:length(true.beta),
                                     which.min(glasso.aic)]
  beta.bic = fit.glasso$coefficients[1+0:length(true.beta),
                                     which.min(glasso.bic)]
  beta = list(beta.aic,beta.bic)

}
else if(method == 'glasso-PLH')
{
  if(missing(df)){
    df = as.integer(nn**(1/3))
  }
  #GX = 1-ecdf(C)(X-1e-8)
  beta.CWY = PO.ipcw(Surv(X, delta), Z,GX)
  fit.bs = PO.it(Surv(X, delta), Z, df = df,
                 beta = beta.CWY,
                 # lam = rep(0,10),
                 GX = GX,
                 glm.maxit = 1, lam.type = "newton")
  fit.glasso = PO.glasso(delta, Z, fit.bs$bs.x,
                         grp,
                         fit.bs$beta, fit.bs$lam, fit.bs$a)
  glasso.refit = PO.bs.ploglik(Surv(X, delta), Z, fit.bs$bs.x,
                               fit.glasso$coefficients[1+1:ncol(dat$Z),],
                               fit.bs$lam, GX
  )
  glasso.ploglik = glasso.refit$loglik*nn


  glasso.paic = 2*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.ploglik
  glasso.pbic = log(nn)*(apply(fit.glasso$coefficients!=0, 2, sum)-df-1) - 2*glasso.ploglik

  beta.paic = fit.glasso$coefficients[1+0:length(true.beta),
                                      which.min(glasso.paic)]
  beta.pbic = fit.glasso$coefficients[1+0:length(true.beta),
                                      which.min(glasso.pbic)]
  beta = list(beta.paic,beta.pbic)
}
}

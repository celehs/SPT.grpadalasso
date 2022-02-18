PO.fit = function (delta, X, Z, C, df  , order = 1,
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

  if(method == 'U-method'){
    beta = PO.ipcw(Surv(X, delta), Z,GX, maxit = control$ipcw.maxit,tol = control$ipcw.tol)
  }
  else if(method == 'B-spline'){
    if(missing(df)){
      df = as.integer(nn**(1/3))
    }
    beta.CWY = PO.ipcw(Surv(X, delta), Z,GX,
                       maxit = control$ipcw.maxit,
                       tol = control$ipcw.tol)
    # Run my algorithm
    fit.bs = PO.it(Surv(X, delta), Z, df = df,
                   beta = beta.CWY,
                   # lam = rep(0,10),
                   GX = GX,
                   glm.maxit = control$it.glm.maxit,
                   lam.type = control$it.lam.type,
                   maxit = control$it.maxit,
                   tol = control$it.tol,
                   max.move = control$it.max.move,
                   order = control$order,
                   )

    # Check estimation of beta
    beta = fit.bs$beta
  }
  else if(method == 'NPMLE')
  {
    beta.CWY = PO.ipcw(Surv(X, delta), Z,GX,
                       maxit = control$ipcw.maxit,
                       tol = control$ipcw.tol)
    tseq = sort(unique(X[delta==1]))
    ht = PO.base.ipcw(tseq,Surv(X, delta), Z,
                      1-ecdf(C)(tseq-1e-8),
                      beta.CWY,
                      maxit = control$base.ipcw.maxit,
                      tol = control$base.ipcw.tol,
                      max.move = control$max.move)
    fit.npmle = PO.NPMLE(Surv(X, delta), Z,
                         beta = beta.CWY,
                         lam = c(0,diff(ht)),
                         maxit = control$NPMLE.maxit,
                         tol = control$NPMLE.tol,
                         glm.maxit = control$NPMLE.glm.maxit,
                         max.move = control$NPMLE.max.move)
    beta = fit.npmle$beta
  }

  res = list()
  res$coefficients = beta
  return(res)
  print(beta)
}

PO(1,1,1,1,method='z')

PO.fit = function (delta, X, Z, C, df  , order = 1,
               method = c('U-method','B-spline','NPMLE'),control)
{
  nn = length(delta)
  tseq = sort(unique(X[delta==1]))
  method <- match.arg(method)
  res = list()
  if (missing(C)){
    KMfit = summary(survfit(Surv(X,delta)~1), time = sort(unique(X)))
    GX = stepfun(KMfit$time, c(1,KMfit$surv))(X)
    Gtseq = stepfun(KMfit$time, c(1,KMfit$surv))(tseq)
  }
  else{
    GX = 1-ecdf(C)(X-1e-8)
    Gtseq = 1-ecdf(C)(tseq-1e-8)
  }

  if(method == 'U-method'){
    beta = PO.ipcw(Surv(X, delta), Z,GX, maxit = control$ipcw.maxit,tol = control$ipcw.tol)


    ht = PO.base.ipcw(tseq,Surv(X, delta), Z,
                      Gtseq,
                      beta)
    ht[is.infinite(ht)] = min(ht[!is.infinite(ht)])
    baseline = stepfun(tseq,c(ht[1],ht))
    res$baseline = baseline
    res$coefficients = beta
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

    res = fit.bs
    names(res)[1] = 'coefficients'
  }
  else if(method == 'NPMLE')
  {

    beta.CWY = PO.ipcw(Surv(X, delta), Z,GX,
                       # maxit = control$ipcw.maxit,
                       # tol = control$ipcw.tol
                       )
    tseq = sort(unique(X[delta==1]))

    ht = PO.base.ipcw(tseq,Surv(X, delta), Z,
                     # 1-ecdf(C)(tseq-1e-8),
                      Gtseq,
                      beta.CWY,
                      maxit = control$base.ipcw.maxit,
                      tol = control$base.ipcw.tol,
                      max.move = control$max.move
                     )

    ht[is.infinite(ht)] = min(ht[!is.infinite(ht)])
    fit.npmle = PO.NPMLE(Surv(X, delta), Z,
                         beta = beta.CWY,
                         lam = c(0,diff(ht)),
                         maxit = control$NPMLE.maxit,
                         tol = control$NPMLE.tol,
                         glm.maxit = control$NPMLE.glm.maxit,
                         max.move = control$NPMLE.max.move
                         )
    res$coefficients = fit.npmle$beta
  }
  return(res)
}

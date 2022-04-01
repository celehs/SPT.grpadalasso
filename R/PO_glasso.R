require(grplasso)

# delta = dat$delta
# Z = dat$Z
# bs.x = fit.bs$bs.x
# beta.mle = fit.bs$beta
# a.mle = fit.bs$a
# lam.mle = fit.bs$lam
# grp = c(1,2,1)
# lambda.min = 1e-4
# nlambda = 100
# min.ratio = 0.95

PO.glasso = function(delta, Z, bs.x, grp,
                     beta.mle, lam.mle, a.mle,
                     ada = TRUE,
                     lambda.list,
                     lambda.min = 1e-4, nlambda = 100)
{
  n = length(delta)
  p = ncol(Z)
  df = ncol(bs.x$Ibs)
  lp.mle = drop(Z %*% beta.mle) + a.mle + drop(bs.x$Ibs %*% lam.mle)
  wgt = c((1+delta)* expit(lp.mle)*expit(-lp.mle),
          drop(delta %*% bs.x$bs)/n/lam.mle^2 )

  X.pseudo = rbind(cbind(1,Z,bs.x$Ibs),
                   cbind(matrix(0,df,p+1), diag(1,df,df)))
  Y.pseudo = c(lp.mle, lam.mle)

  index = c(NA, grp, rep(NA, df))

  scale.factor = rep(1, length(index))
  if(ada)
  {
    for(igrp in unique(index))
    {
      if(is.na(igrp))
        next
      igrp.pos = which(igrp == index)
      grp.l2 = sqrt(sum(beta.mle[igrp.pos-1]^2))
      X.pseudo[,igrp.pos] = X.pseudo[,igrp.pos]*grp.l2
      scale.factor[igrp.pos] = grp.l2
    }
  }

  if(missing(lambda.list))
  {
    lambda.max = lambdamax(X.pseudo, y = Y.pseudo,
                           weights = wgt,
                           index = index,
                       model = LinReg(), center = FALSE, standardize=FALSE)
    lambda.list = (seq((lambda.max),
                          (lambda.min),
                          length.out = nlambda))
  }
  # print(scale.factor)
  sink("NUL")
  fit = grplasso(X.pseudo, y = Y.pseudo,
                 weights = wgt,
                 index = index,
                 model = LinReg(), center = FALSE, standardize=FALSE,
                 lambda = lambda.list)
  sink()
  fit$LSA = apply(wgt*(Y.pseudo - X.pseudo %*% fit$coefficients)^2,2,mean)
  fit$coefficients = fit$coefficients * scale.factor
  fit$scale.factor = scale.factor
  return(fit)
}

# PO.glasso.cv = function(surv, Z, C, bs.x, grp,
#                                beta.mle, lam.mle, a.mle,
#                                ada = TRUE, lambda.list,
#                                nfold = 5, tseq,
#                             measure = c("ploglik","loglik",
#                                         "APE", "APERefit"))
# {
#   n = length(surv)
#   foldid = rep_len(1:nfold, n)[sample(1:n)]
#
#   if(missing(tseq))
#   {
#     tseq = sort(unique(surv[surv[,2]==1,1]))
#   }
#
#   meas.cv = matrix(NA, nfold, length(lambda.list))
#   for(ifold in 1:nfold)
#   {
#     train.pos = which(foldid != ifold)
#     test.pos = which(foldid == ifold)
#
#     bs.fold = PO.it(surv[train.pos], Z[train.pos,],
#                    bs.knots = bs.x$knots,
#                    beta = beta.mle, lam = lam.mle, a = a.mle,
#                    GX = GX,
#                    glm.maxit = 1, lam.type = "newton")
#     glasso.fold = PO.glasso(surv[train.pos,2], Z[train.pos,], bs.fold$bs.x,
#                             grp,
#                             bs.fold$beta, bs.fold$lam, bs.fold$a,
#                             ada = ada,
#                             lambda.list = lambda.list)
#     bs.test = bs.x
#     bs.test$bs = bs.test$bs[test.pos,]
#     bs.test$Ibs = bs.test$Ibs[test.pos,]
#
#     if(measure[1] == "loglik")
#     {
#       meas.cv[ifold, ] = PO.bs.loglik(surv[test.pos,2], Z[test.pos,], bs.test,
#                    glasso.fold$coefficients)
#     }else if(measure[1] == "ploglik")
#     {
#       meas.cv[ifold, ] = PO.bs.ploglik(surv[test.pos], Z[test.pos,],
#                                        bs.test,
#                                        glasso.fold$coefficients[1+1:ncol(Z),],
#                                        bs.fold$lam, GX)$loglik
#     }else if(measure[1] == "APE")
#     {
#       for(i in 1:length(lambda.list))
#       {
#         lp.bs = (drop(Z[test.pos,] %*% fit.adaglasso$coefficients[1+1:ncol(Z),i]) +
#                    fit.adaglasso$coefficients[1,i])
#         htseq = Ibs.1.pred(bs.x$knots, tseq,
#                            glasso.fold$coefficients[-(1+0:ncol(Z)),i])
#         valid.tmp = expit(outer(lp.bs,htseq,'+'))
#
#         DT = pred.DT(C[test.pos], tseq, valid.tmp)
#         thres.list = quantile(DT$piC, probs = seq(0,1,0.01))
#         Xhat = pred.Xhat(C[test.pos], DT, thres.list)
#         APE = pred.APE(surv[test.pos,1], Xhat, thres.list)
#         meas.cv[ifold, i] = APE$best
#       }
#     }else if(measure[1] == "APERefit")
#     {
#       refit = PO.bs.ploglik(surv[train.pos], Z[train.pos,],
#                             bs.fold$bs.x,
#                             glasso.fold$coefficients[1+1:ncol(Z),],
#                             bs.fold$lam, GX)
#       for(i in 1:length(lambda.list))
#       {
#         lp.bs = (drop(Z[test.pos,] %*% fit.adaglasso$coefficients[1+1:ncol(Z),i]) +
#                    refit$a[i])
#         htseq = Ibs.1.pred(bs.x$knots, tseq,
#                            refit$lam[,i])
#         valid.tmp = expit(outer(lp.bs,htseq,'+'))
#
#         DT = pred.DT(C[test.pos], tseq, valid.tmp)
#         thres.list = quantile(DT$piC, probs = seq(0,1,0.01))
#         Xhat = pred.Xhat(C[test.pos], DT, thres.list)
#         APE = pred.APE(surv[test.pos,1], Xhat, thres.list)
#         meas.cv[ifold, i] = APE$best
#       }
#     }
#
#   }
#
#   return(meas.cv)
# }

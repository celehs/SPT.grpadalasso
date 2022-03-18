
PO.glasso.control = function (inv.m = exp,  #PO.sim
                       ipcw.maxit = 1000, ipcw.tol = 1e-6,   #PO.ipcw
                       #PO.it
                       it.maxit = 1000, it.tol = 1e-6,it.glm.maxit = 1, it.max.move = 1, order = 1,it.lam.type = c("newton","coord-joint", "coord-loop","breslow"),
                       #glasso
                       ##PO.glasso
                       ada = TRUE, lambda.min = 1e-4, nlambda = 100,
                       ##PO.bs.ploglik
                       bs.ploglik.maxit = 1000, bs.ploglik.tol = 1e-6,bs.ploglik.glm.maxit = 1, bs.ploglik.max.move = 1, bs.ploglik.lam.type = c("newton","coord-joint", "coord-loop","breslow"),
                       #origin
                       iter.max = 20, outer.max = 10, timefix = TRUE)
{
  if (min(ipcw.maxit,it.maxit,bs.ploglik.maxit) < 0)
    stop("Invalid value for iterations")
  if (min(ipcw.tol,it.tol,bs.ploglik.tol) < 0)
    stop("Invalid value for tolerance")
  if (min(it.glm.maxit,bs.ploglik.glm.maxit) < 0)
    stop("Invalid value for glm.maxit")
  if (min(it.max.move,bs.ploglik.max.move) < 0)
    stop("Invalid value for max.move")
  it.lam.type <- match.arg(it.lam.type)
  bs.ploglik.lam.type <- match.arg(bs.ploglik.lam.type)

  list(inv.m = inv.m, ipcw.maxit = ipcw.maxit, ipcw.tol = ipcw.tol,
       it.maxit = it.maxit, it.tol = it.tol,it.glm.maxit = it.glm.maxit, it.max.move = it.max.move, order = order,it.lam.type = it.lam.type,
       ada = ada, lambda.min = lambda.min, nlambda = nlambda,
       bs.ploglik.maxit = bs.ploglik.maxit, bs.ploglik.tol = bs.ploglik.tol,bs.ploglik.glm.maxit = bs.ploglik.glm.maxit, bs.ploglik.max.move = bs.ploglik.max.move,
       bs.ploglik.lam.type = bs.ploglik.lam.type,
       iter.max = iter.max, outer.max = outer.max, timefix = timefix

  )
}


PO.control = function (inv.m = exp,  #PO.sim
                       ipcw.maxit = 1000, ipcw.tol = 1e-6,   #PO.ipcw
                       #PO.it
                       it.maxit = 1000, it.tol = 1e-6,it.glm.maxit = 1, it.max.move = 1, order = 1,it.lam.type = c("newton","coord-joint", "coord-loop","breslow"),
                       #NPMLE
                       ##PO.base.ipcw
                       base.ipcw.maxit = 1000, base.ipcw.tol = 1e-8, base.ipcw.max.move = 1,
                       ##PO.NPMLE
                       NPMLE.maxit = 1000, NPMLE.tol = 1e-6,NPMLE.glm.maxit = 1, NPMLE.max.move = 1
                        )
{
  if (min(ipcw.maxit,it.maxit,base.ipcw.maxit,NPMLE.maxit) < 0)
    stop("Invalid value for iterations")
  if (min(ipcw.tol,it.tol,base.ipcw.tol,NPMLE.tol) < 0)
    stop("Invalid value for tolerance")
  if (min(it.glm.maxit,NPMLE.glm.maxit) < 0)
    stop("Invalid value for glm.maxit")
  if (min(it.max.move,base.ipcw.max.move,NPMLE.max.move) < 0)
    stop("Invalid value for max.move")
  it.lam.type <- match.arg(it.lam.type)

  list(inv.m = inv.m, ipcw.maxit = ipcw.maxit, ipcw.tol = ipcw.tol,
       it.maxit = it.maxit, it.tol = it.tol,it.glm.maxit = it.glm.maxit, it.max.move = it.max.move, order = order,it.lam.type = it.lam.type,
       base.ipcw.maxit = base.ipcw.maxit, base.ipcw.tol = base.ipcw.tol, base.ipcw.max.move = base.ipcw.max.move,
       NPMLE.maxit = NPMLE.maxit, NPMLE.tol = NPMLE.tol,NPMLE.glm.maxit = NPMLE.glm.maxit, NPMLE.max.move = NPMLE.max.move
       )
}

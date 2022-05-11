require(survival)

# dat = PO.sim(10000, c(0,1,-0.5))
# surv = Surv(dat$X, dat$delta)
# Z = dat$Z
# df = 10
# maxit = 1000
# tol = 1e-6
# glm.maxit = 1
# beta = rep(0, p)
# lam = rep(0, length(sort(unique(surv[surv[,2]==1,1]))))

# surv = Surv(SX, Delta)
# beta = bb.init
# maxit = 1000
# tol = 1e-6
# glm.maxit = 1
# Z = Z.old
# max.move = 1
# a = log(tol)
# lam = pmax(0,c(0,diff(ht)))

PO.NPMLE = function(surv, Z,
                 beta, lam,  maxit = 1000, tol = 1e-6,
                 glm.maxit = 1, max.move = 1, a = log(tol))
{

  bf.t1 = surv[,1] < min(surv[surv[,2]==1,1])
  surv = surv[!bf.t1]
  Z = Z[!bf.t1,]

  n = length(surv)
  p = ncol(Z)

  # Calculate the baseline
  # Computation can be improved by order the observation times
  X.order = order(surv[,1])
  surv = surv[X.order]
  Z = Z[X.order,]
  tk = surv[surv[,2]==1,1]
  dup.tk = duplicated(tk)
  lam.num = diff(c((1:length(tk))[!dup.tk], length(tk)+1)) # ties at tk
  tk = tk[!dup.tk]
  tk.pos = match(tk, surv[,1])
  X.to.tk = rep( 0:length(tk),diff(c(0,tk.pos-1,n)))+1

  # Initial value
  if(missing(beta))
  {
    beta = rep(0, p)
  }
  if(missing(lam))
  {
    lam = rep(0, length(tk))
  }
  lp = drop(Z %*% beta)
  mX = c(0,cumsum(lam))[X.to.tk]

  # Pseudo data
  d.pseudo = c(rep(0:1, c(n,sum(surv[,2]))))
  pseudo.pos = c(1:n , which(surv[,2]==1))
  Z.pseudo = Z[pseudo.pos,]

  # Iterative algorithm
  # lam.num = drop(surv[,2] %*% X.eq.t)
  iter = 0

  cont.lam = TRUE
  while (cont.lam)
  {
    a.new = coef(glm(d.pseudo ~ 1, family = binomial,
                 offset = (mX+lp-a)[pseudo.pos]))

    lam.new =  lam.num/(rev(cumsum(rev((1+surv[,2])*expit(lp+mX)-surv[,2])))[tk.pos])
    # neg.lam = lam.new <0
    # lam.new[neg.lam] = lam[neg.lam]/2
    update = lam.new - lam
    rescale =  min(1, max.move/sqrt(sum(update^2+(a.new -a)^2)))
    lam.new = lam + update *rescale
    a.new = a + (a.new-a) *rescale
    iter = iter + 1
    if(iter > maxit)
    {
      stop("Fail to converge at max iteration.")
    }
    cont.lam =max(abs(expit(cumsum(c(a,lam)))
                   - expit(cumsum(c(a.new,lam.new))))) > tol

    lam = lam.new
    lp = lp - a + a.new
    a = a.new
    mX = c(0,cumsum(lam))[X.to.tk]
  }

  cont.beta = TRUE
  while(cont.beta)
  {
    a.new = coef(glm(d.pseudo ~ 1, family = binomial,
                     offset = (mX+lp-a)[pseudo.pos]))
    lam.new =  lam.num/(rev(cumsum(rev((1+surv[,2])*expit(lp+mX)-surv[,2])))[tk.pos])
    # neg.lam = lam.new <0
    # lam.new[neg.lam] = lam[neg.lam]/2
    update = lam.new - lam
    rescale =  min(1, max.move/sqrt(sum(update^2+(a.new -a)^2)))
    lam.new = lam + update *rescale
    a.new = a + (a.new-a) *rescale

    mX = c(0,cumsum(lam.new))[X.to.tk]
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
    update = coef(tmp.fit)-c(a.new,beta)
    rescale =  min(1, max.move/sqrt(sum(update^2)))
    iter = iter + tmp.fit$iter
    cont.beta = max(abs(update[-1]),
                    abs(expit(cumsum(c(a,lam)))
                        - expit(cumsum(c(a.new+update[1],lam.new))))) > tol
    beta = beta + update[-1]*rescale
    a = a.new+update[1]
    lp = a+drop(Z %*% beta)
    lam = lam.new
    if(iter > maxit)
    {
      stop("Fail to converge at max iteration.")
    }
  }

  return(list(beta = beta,  lam = lam, a = a,
              tk = tk, iter = iter))

}

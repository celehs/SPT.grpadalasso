#' NPMLE for Logistic-CoxPH Cure-Rate Model
#' @usage PO(formula, data, C, df, weights, subset, init,control,
#' singular.ok = TRUE,model = FALSE,x = FALSE, y = TRUE, tt,
#' method = c('U-method','B-spline','NPMLE','glasso','glasso-PLH'),...)
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the \code{Surv.cure} function.
#' @param data a \code{data.frame} in which to interpret the variables named in
#'  the \code{formula}, or in the \code{subset} and the \code{weights} argument.
#' @param C Censoring time
#' @param df degree of freedom
#' @param weights vector of case weights, see the note below. For a thorough
#'  discussion of these see the book by Therneau and Grambsch.
#' @param subset expression indicating which subset of the rows of data should
#'  be used in the fit. All observations are included by default.
#' @param init vector of initial values of the iteration. Default initial value
#'   is zero for all variables.
#' @param control Object of class \code{PO.control} specifying iteration
#'   limit and other control options. Default is \code{PO.control(...)}.
#' @param singular.ok   logical value indicating how to handle collinearity in
#'  the model matrix. If \code{TRUE}, the program will automatically skip over
#'  columns of the X matrix that are linear combinations of earlier columns.
#'  In this case the coefficients for such columns will be \code{NA}, and the
#'  variance matrix will contain zeros. For ancillary calculations, such as the
#'  linear predictor, the missing coefficients are treated as zeros.
#' @param model logical value: if \code{TRUE}, the model frame is returned in
#'  component model.
#' @param x logical value: if \code{TRUE}, the x matrix is returned in
#'  component \code{x}.
#' @param y logical value: if \code{TRUE}, the response vector is returned in component \code{y}.
#' @param tt optional list of time-transform functions.
#' @param method a character string specifying the method in \code{c('U-method',
#' 'B-spline','NPMLE')} for estimation. The default method is the \code{U-method}.
#' @param ... other parameters passed to \code{PO.control}.
#' @author Jue Hou
#' @examples
#'  # A simulated data set
#' require(survival)
#' data('sim_PO_data')
#' res = PO(Surv(X, delta) ~ Z[,1]+ Z[,2]+ Z[,3],data = sim_PO_data)

#' # Or you may generate another one
#' set.seed(1)
#' df = 10
#' nn = 1000
#' beta = c(0.5,0,-0.5, rep(0,10))
#' sim_PO_data = PO.sim(nn, beta,
#'                      C.gen = function(n)
#'                        5+rbinom(n,1,0.5)*runif(n, -5, 0))

#' # Fit PO model
#' res = PO(Surv(X, delta) ~ Z[,1]+ Z[,2]+ Z[,3]+ Z[,4]+ Z[,5]+ Z[,6],data = sim_PO_data)
#' @export
PO = function (formula, data, C, df, weights, subset, init,control,
           singular.ok = TRUE,model = FALSE,
           x = FALSE, y = TRUE, tt,
           method = c('U-method','B-spline','NPMLE','glasso','glasso-PLH'),...)
{
  Call <- match.call()
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(PO.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(indx == 0L))
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx ==
                                                                  0L]), domain = NA)
  }
  method <- match.arg(method)
  if (missing(control))
  {
    if(method %in% c('U-method','B-spline','NPMLE'))
    {
      control <- PO.control(...)
    }
    else if(method %in% c('glasso','glasso-PLH'))
    {
      control <- PO.glasso.control(...)
    }
  }
  if (missing(formula))
    stop("a formula argument is required")
  ss <- c("cluster", "offset")
  if (is.list(formula))
    Terms <- if (missing(data))
      terms(formula[[1]], specials = ss)
  else terms(formula[[1]], specials = ss, data = data)
  else Terms <- if (missing(data))
    terms(formula, specials = ss)
  else terms(formula, specials = ss, data = data)
  tcl <- attr(Terms, "specials")$cluster
  if (length(tcl) > 1)
    stop("a formula cannot have multiple cluster terms")
  if (length(tcl) > 0) {
    factors <- attr(Terms, "factors")
    if (any(factors[tcl, ] > 1))
      stop("cluster() cannot be in an interaction")
    if (attr(Terms, "response") == 0)
      stop("formula must have a Surv response")
    temp <- attr(Terms, "term.labels")
    oo <- attr(Terms, "specials")$offset
    if (!is.null(oo)) {
      ooterm <- rownames(factors)[oo]
      if (oo < tcl)
        temp <- c(ooterm, temp)
      else temp <- c(temp, ooterm)
    }
    if (is.null(Call$cluster))
      Call$cluster <- attr(Terms, "variables")[[1 + tcl]][[2]]
    else warning("cluster appears both in a formula and as an argument, formula term ignored")
    if (is.list(formula))
      formula[[1]][[3]] <- reformulate(temp[1 - tcl])[[2]]
    else formula[[3]] <- reformulate(temp[1 - tcl])[[2]]
    Call$formula <- formula
  }
  indx <- match(c("formula", "data", "weights", "subset"), names(Call),
                nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")
  tform <- Call[c(1, indx)]
  tform[[1L]] <- quote(stats::model.frame)
  if (is.list(formula)) {
    multiform <- TRUE
    dformula <- formula[[1]]
    tlab <- unlist(lapply(covlist$rhs, function(x) attr(terms.formula(x$formula),
                                                        "term.labels")))
    tlab <- c(attr(terms.formula(dformula), "term.labels"),
              tlab)
    newform <- reformulate(tlab, dformula[[2]])
    environment(newform) <- environment(dformula)
    formula <- newform
  }
  else {
    multiform <- FALSE
    covlist <- NULL
    dformula <- formula
  }
  special <- c()
  tform$formula <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data = data)
  if (!is.null(attr(tform$formula, "specials")$tt)) {
    coxenv <- new.env(parent = environment(formula))
    assign("tt", function(x) x, envir = coxenv)
    environment(tform$formula) <- coxenv
  }
  mf <- eval(tform, parent.frame())
  Terms <- terms(mf)
  n <- nrow(mf)
  Y <- model.response(mf)
  isSurv2 <- inherits(Y, "Surv2")
  if (isSurv2) {
    if (!is.null(attr(Terms, "specials")$cluster))
      stop("cluster() cannot appear in the model statement")
    new <- Surv2data(mf)
    mf <- new$mf
    Y <- new$y
    n <- nrow(mf)
  }
  else {
    if (!is.Surv(Y))
      stop("Response must be a survival object")
  }
  if (n == 0)
    stop("No (non-missing) observations")
  type <- attr(Y, "type")
  multi <- FALSE
  if (type == "mright" || type == "mcounting")
    multi <- TRUE
  else if (type != "right" && type != "counting")
    stop(paste("Cox model doesn't support \"", type, "\" survival data",
               sep = ""))
  data.n <- nrow(Y)
  if (!multi && multiform)
    stop("formula is a list but the response is not multi-state")
  if (multi && length(attr(Terms, "specials")$frailty) > 0)
    stop("multi-state models do not currently support frailty terms")
  if (multi && length(attr(Terms, "specials")$pspline) > 0)
    stop("multi-state models do not currently support pspline terms")
  if (multi && length(attr(Terms, "specials")$ridge) > 0)
    stop("multi-state models do not currently support ridge penalties")
  strats <- attr(Terms, "specials")$strata
  hasinteractions <- FALSE
  dropterms <- NULL
  if (length(strats)) {
    stemp <- untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) == 1)
      strata.keep <- mf[[stemp$vars]]
    else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
    istrat <- as.integer(strata.keep)
    for (i in stemp$vars) {
      if (any(attr(Terms, "order")[attr(Terms, "factors")[i,
      ] > 0] > 1))
        hasinteractions <- TRUE
    }
    if (!hasinteractions)
      dropterms <- stemp$terms
  }
  else istrat <- NULL
  if (hasinteractions && multi)
    stop("multi-state coxph does not support strata*covariate interactions")
  timetrans <- attr(Terms, "specials")$tt
  if (missing(tt))
    tt <- NULL
  if (length(timetrans)) {
    if (multi || isSurv2)
      stop("the tt() transform is not implemented for multi-state or Surv2 models")
    timetrans <- untangle.specials(Terms, "tt")
    ntrans <- length(timetrans$terms)
    if (is.null(tt)) {
      tt <- function(x, time, riskset, weights) {
        obrien <- function(x) {
          r <- rank(x)
          (r - 0.5)/(0.5 + length(r) - r)
        }
        unlist(tapply(x, riskset, obrien))
      }
    }
    if (is.function(tt))
      tt <- list(tt)
    if (is.list(tt)) {
      if (any(!sapply(tt, is.function)))
        stop("The tt argument must contain function or list of functions")
      if (length(tt) != ntrans) {
        if (length(tt) == 1) {
          temp <- vector("list", ntrans)
          for (i in 1:ntrans) temp[[i]] <- tt[[1]]
          tt <- temp
        }
        else stop("Wrong length for tt argument")
      }
    }
    else stop("The tt argument must contain a function or list of functions")
    if (ncol(Y) == 2) {
      if (length(strats) == 0) {
        sorted <- order(-Y[, 1], Y[, 2])
        newstrat <- rep.int(0L, nrow(Y))
        newstrat[1] <- 1L
      }
      else {
        sorted <- order(istrat, -Y[, 1], Y[, 2])
        newstrat <- as.integer(c(1, 1 * (diff(istrat[sorted]) !=
                                           0)))
      }
      if (storage.mode(Y) != "double")
        storage.mode(Y) <- "double"
    }
    else {
      if (length(strats) == 0) {
        sort.end <- order(-Y[, 2], Y[, 3])
        sort.start <- order(-Y[, 1])
        newstrat <- c(1L, rep(0, nrow(Y) - 1))
      }
      else {
        sort.end <- order(istrat, -Y[, 2], Y[, 3])
        sort.start <- order(istrat, -Y[, 1])
        newstrat <- c(1L, as.integer(diff(istrat[sort.end]) !=
                                       0))
      }
      if (storage.mode(Y) != "double")
        storage.mode(Y) <- "double"
    }
    type <- "right"
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights)))
      stop("weights must be finite")
    tcall <- attr(Terms, "variables")[timetrans$terms +
                                        2]
    pvars <- attr(Terms, "predvars")
    pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
    for (i in 1:ntrans) {
      newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[, 1],
                         istrat, weights)
      mf[[timetrans$var[i]]] <- newtt
      nclass <- class(newtt)
      if (any(nclass %in% pmethod)) {
        dummy <- as.call(list(as.name(class(newtt)[1]),
                              tcall[[i]][[2]]))
        ptemp <- makepredictcall(newtt, dummy)
        pvars[[timetrans$terms[i] + 2]] <- ptemp
      }
    }
    attr(Terms, "predvars") <- pvars
  }
  xlevels <- .getXlevels(Terms, mf)
  cluster <- model.extract(mf, "cluster")
  weights <- model.weights(mf)
  has.cluster <- !(missing(cluster) || length(cluster) ==
                     0)
  has.rwt <- (!is.null(weights) && any(weights != floor(weights)))
  ncluster <- 0
  contrast.arg <- NULL
  attr(Terms, "intercept") <- 1
  if (length(dropterms)) {
    Terms2 <- Terms[-dropterms]
    X <- model.matrix(Terms2, mf, constrasts.arg = contrast.arg)
    temp <- attr(X, "assign")
    shift <- sort(dropterms)
    for (i in seq(along.with = shift)) temp <- temp + 1 *
      (shift[i] <= temp)
    attr(X, "assign") <- temp
  }
  else X <- model.matrix(Terms, mf, contrasts.arg = contrast.arg)
  Xatt <- attributes(X)
  if (hasinteractions)
    adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  else adrop <- 0
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset == 0))
    offset <- rep(0, nrow(mf))
  else if (any(!is.finite(offset) | !is.finite(exp(offset))))
    stop("offsets must lead to a finite risk score")
  weights <- model.weights(mf)
  if (!is.null(weights) && any(!is.finite(weights)))
    stop("weights must be finite")
  assign <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  if (sum(Y[, ncol(Y)]) == 0) {
    ncoef <- ncol(X)
    ctemp <- rep(NA, ncoef)
    names(ctemp) <- colnames(X)
    concordance = c(concordant = 0, discordant = 0, tied.x = 0,
                    tied.y = 0, tied.xy = 0, concordance = NA, std = NA,
                    timefix = FALSE)
    rval <- list(coefficients = ctemp, var = matrix(0, ncoef,
                                                    ncoef), loglik = c(0, 0), score = 0, iter = 0, linear.predictors = offset,
                 residuals = rep(0, data.n), means = colMeans(X),
                 #method = method,
                 n = data.n, nevent = 0, terms = Terms,
                 assign = assign, concordance = concordance, wald.test = 0,
                 y = Y, call = Call)
    class(rval) <- "coxph"
    return(rval)
  }
  if (!all(is.finite(X)))
    stop("data contains an infinite predictor")
   if (missing(init))
     init <- NULL
  else {
    if (length(init) != ncol(X))
      stop("wrong length for init argument")
    temp <- X %*% init - sum(colMeans(X) * init) + offset
    if (any(exp(temp) > .Machine$double.xmax) || all(exp(temp) ==
                                                     0))
      stop("initial values lead to overflow or underflow of the exp function")
  }
  nn = length(Y[,2])
  if(missing(df)){
    df = as.integer(nn**(1/3))
  }

  if(method %in% c('U-method','B-spline','NPMLE'))
  {

    if (missing(C)){
      fit = PO.fit(delta = Y[,2], X = Y[,1], Z = X, df = df  , order = 1,
                   method = method,control = control)
    }
    else{
      fit = PO.fit(delta = Y[,2], X = Y[,1], Z = X, C = C, df = df  , order = 1,
                   method = method,control = control)
    }
  }
  else if(method %in% c('glasso','glasso-PLH'))
  {
    if (missing(C)){
      fit = PO.glasso.fit(delta = Y[,2], X = Y[,1], Z = X, df = df  , order = 1,
                   method = method,control = control)
    }
    else{
      fit = PO.glasso.fit(delta = Y[,2], X = Y[,1], Z = X, C = C, df = df  , order = 1,
                   method = method,control = control)
    }
  }
  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "PO"
  }
  else {
    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable",
                   paste(vars, collapse = " "))
      if (!singular.ok)
        stop(msg)
    }
    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)])
    fit$terms <- Terms
    class(fit) <- fit$class
    fit$class <- NULL

    if (model) {
      if (length(timetrans)) {
        stop("'model=TRUE' not supported for models with tt terms")
      }
      fit$model <- mf
    }
    if (x) {
      fit$x <- X
    }
    if (y)
      fit$y <- Y
  }
  fit$formula <- formula(Terms)
  fit$call <- Call
  fit
}

# GX = 1-ecdf(sim_PO_data$C)(sim_PO_data$X-1e-8)
# beta.CWY = PO.ipcw(Surv(sim_PO_data$X, sim_PO_data$delta), sim_PO_data$Z,GX)

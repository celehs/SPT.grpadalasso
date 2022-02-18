\name{PO.control}
\alias{PO.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancillary arguments for controling PO fits
}
\description{
Assign the ancillary arguments for controling cureph fits
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PO.control(inv.m = exp, ipcw.maxit = 1000, ipcw.tol = 1e-6,
                    it.maxit = 1000, it.tol = 1e-6,it.glm.maxit = 1,
                    it.max.move = 1, order = 1,
                    it.lam.type = c("newton","coord-joint", "coord-loop","breslow"),
                    base.ipcw.maxit = 1000, base.ipcw.tol = 1e-8,
                    base.ipcw.max.move = 1,NPMLE.maxit = 1000, NPMLE.tol = 1e-6,
                    NPMLE.glm.maxit = 1, NPMLE.max.move = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inv.m}{ Sample size of the data}
  \item{ipcw.maxit}{Maximum number of iterations to attempt for convergence}
  \item{ipcw.tol}{Tolerance criteria}
  \item{it.maxit}{Maximum number of iterations to attempt for convergence}
  \item{it.tol}{Tolerance criteria}
  \item{it.glm.maxit}{Maximum number of iterations to attempt for convergence}
  \item{it.max.move}{}
  \item{order}{}
  \item{it.lam.type}{a character string specifying the method for solving.}
  \item{base.ipcw.maxit}{Maximum number of iterations to attempt for convergence}
  \item{base.ipcw.tol}{Tolerance criteria}
  \item{base.ipcw.max.move}{}
  \item{NPMLE.maxit}{Maximum number of iterations to attempt for convergence for NPMLE method}
  \item{NPMLE.tol}{Tolerance criteria for NPMLE method}
  \item{NPMLE.glm.maxit}{Maximum number of iterations to attempt for convergence}
  \item{NPMLE.max.move}{}
}
\details{
\code{line.search} and \code{init.step} only matter when \code{method = "BFGS-Newton"}.
}
\value{
a list containing the values of each of the above constants
}
\author{
Jue (Marquis) Hou
}

\seealso{
\code{\link{cureph}}
}

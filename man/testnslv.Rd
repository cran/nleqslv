\name{testnslv}
\alias{testnslv}
\encoding{UTF-8}
\title{Test different methods for solving with \code{nleqslv}}
\description{
The function tests different methods and global strategies for solving a system of nonlinear equations with \code{nleqslv}
}
\usage{
testnslv(x, fn, jac=NULL, ...,
          method = c("Newton", "Broyden"),
          global = c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"),
          Nrep=0L, title=NULL
        )
}
\arguments{
 \item{x}{A numeric vector with an initial guess of the root.}
 \item{fn}{A function of \code{x} returning the function values.}
 \item{jac}{A function to return the Jacobian for the \code{fn} function.
            For a vector valued function \code{fn} the Jacobian must be a numeric
            matrix of the correct dimensions.
            For a scalar valued function \code{fn} the \code{jac} function may return a scalar.
            If not supplied numerical derivatives will be used.}
 \item{\dots}{Further arguments to be passed to \code{fn} and \code{jac} and \code{\link{nleqslv}}.}
 \item{method}{The methods to use for finding a solution.}
 \item{global}{The global strategies to test. The argument may consist of several possibly abbreviated items.}
 \item{Nrep}{Number of repetitions to apply. Default is no repetitions.}
 \item{title}{a description of this experiment.}
}
\details{The function solves the function \code{fn} with \code{\link{nleqslv}} for the specified methods and global strategies.
When argument \code{Nrep} has been set to a number greater than or equal to 1,
repetitions of the solving process are performed and the used CPU time in seconds is recorded.

If checking a user supplied jacobian is enabled, then \code{testnslv} will stop immediately when a possibly
incorrect jacobian is detected.
}
\value{
\code{testnslv} returns an object of class \code{"test.nleqslv"} which is a list containing the following elements
\describe{
    \item{\code{call}}{the matched call}
    \item{\code{out}}{ a dataframe containing the results with the following columns
        \describe{
          \item{\code{Method}}{method used.}
          \item{\code{Global}}{global strategy used.}
          \item{\code{termcd}}{termination code of \code{nleqslv}.}
          \item{\code{Fcnt}}{number of function evaluations used by the method and global strategy.
                             This excludes function evaluations made when computing a numerical Jacobian.}
          \item{\code{Jcnt}}{number of Jacobian evaluations.}
          \item{\code{Iter}}{number of outer iterations used by the algorithm.}
          \item{\code{Message}}{a string describing the termination code in an abbreviated form.}
          \item{\code{Fnorm}}{square of the euclidean norm of the vector of function results divided by 2.}
          \item{\code{cpusecs}}{CPU seconds used by the requested number of repetitions (only present when
          argument \code{Nrep} is not 0).}
        }
    }
    \item{\code{title}}{the description if specified}
}
The abbreviated strings are
\describe{
\item{\code{Fcrit}}{Convergence of function values has been achieved.}
\item{\code{Xcrit}}{This means that the relative distance between two
consecutive x-values is smaller than \code{xtol}.}
\item{\code{Stalled}}{The algorithm cannot find an acceptable new point.}
\item{\code{Maxiter}}{Iteration limit \code{maxit} exceeded.}
\item{\code{Illcond}}{Jacobian is too ill-conditioned.}
\item{\code{Singular}}{Jacobian is singular.}
\item{\code{BadJac}}{Jacobian is unusable.}
\item{\code{ERROR}}{\code{nleqslv} stopped because of a fatal error.}
}
}
\section{Warning}{
Any \code{nleqslv} error message will be displayed immediately and
an error for the particular combination of method and global strategy will be recorded in the final dataframe.
}
\examples{
dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}
xstart <- c(0.5,0.5)
fstart <- dslnex(xstart)
testnslv(xstart,dslnex)
# this will encounter an error
xstart <- c(2.0,0.5)
fstart <- dslnex(xstart)
testnslv(xstart,dslnex)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonlinear}
\keyword{optimize}

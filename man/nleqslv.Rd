\name{nleqslv}
\alias{nleqslv}
\encoding{UTF-8}
\title{Solving systems of non linear equations with Broyden or Newton}
\description{
The function solves a system of non linear equations with either a Broyden or a
full Newton method. It provides linesearch and trust region global strategies
for difficult systems.
}
\usage{
nleqslv(x, fn, jac=NULL,...,
               method = c("Broyden", "Newton"),
               global = c("dbldog", "pwldog", "qline", "gline"),
               xscalm = c("fixed","auto"),
               control = list()
	)
}
\arguments{
 \item{x}{A list or numeric vector of starting estimates.}
 \item{fn}{A function of \code{x} returning the function values.}
 \item{jac}{A function to return the Jacobian for the \code{fn} function.
If not supplied numerical derivatives will be used.}
 \item{\dots}{Further arguments to be passed to \code{fn} and \code{jac}.}
 \item{method}{The method used for finding a solution. See \sQuote{Details}.}
 \item{global}{The global strategy to apply. See \sQuote{Details}.}
 \item{xscalm}{The type of x scaling to use. See \sQuote{Details}.}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
}
\details{
The algorithms implemented in \code{nleqslv} are based on Dennis and Schnabel (1996).

Method \code{Broyden} starts with a computed Jacobian of the function and
then updates this Jacobian after each successful iteration using the so-called Broyden update.
This method often shows super linear convergence towards a solution.
When \code{nleqslv} determines that it cannot
continue with the current Broyden matrix it will compute a new Jacobian.

Method \code{Newton} calculates a Jacobian of the function \code{fn} at each iteration.
Close to a solution this method will usually show quadratic convergence.

Both methods apply a so-called (backtracking) global strategy to find
a better (more acceptable) iterate.  The function criterion used by the algorithm is the
sum of squares of the function values and \dQuote{acceptable} means sufficient decrease
of the current function criterion value compared to that of the previous iteration.
A comprehensive discussion of these issues can be found in Dennis and Schnabel (1996).
Both methods apply a QR-decomposition to the Jacobian as implemented in Lapack.
The Broyden method applies a rank 1 update to the Jacobian at the end of each iteration
and uses the methods implemented in Reichel and Gragg (1990).

The four global strategies that may be specified in the \code{global} argument,  are
\describe{
        \item{\code{qline}}{a quadratic linesearch}
        \item{\code{gline}}{a geometric linesearch}
        \item{\code{dbldog}}{a trust region method using the double dogleg method as described in
                             Dennis and Schnabel (1996)}
        \item{\code{pwldog}}{a trust region method using the Powell dogleg method
                             as developed by Powell (1970).}
}
The double dogleg method is the default global strategy employed by this package.

The parameters \code{x} may be scaled during the search for a zero of \code{fn}.
The  \code{xscalm} argument provides two possibilities for scaling
\describe{
\item{\code{fixed}}{the scaling factors are set to the values supplied in
the \code{control} argument and remain unchanged during the iterations.
The scaling factor of any element of \code{x}  should be set
to the inverse of the typical value of that element of \code{x},
ensuring that all elements of \code{x} are approximately equal in size.}

\item{\code{auto}}{the scaling factors
are calculated from the euclidean norms of the columns of the jacobian matrix.
When a new Jacobian is computed, the scaling values will be set to the
euclidean norm of the corresponding column if that is larger than the current scaling value.
Thus the scaling values will not decrease during the iteration.
This is the method described in \enc{Moré}{More} (1978).
Usually manual scaling is preferable.}
}

The algorithm cannot cope with a ill-conditioned or singular Jacobian, in contrast
to previous versions of this package. When a Jacobian becomes singular or very ill-conditioned
the algorithm will return an error.  A Jacobian is considered to be very ill-conditioned
when the estimated inverse condition is less than  or equal to the machine precision.

When the function to be solved returns non-finite values for a parameter vector \code{x}
and the algorithm is \emph{not} evaluating a numerical jacobian, then the non-finite
value will be replaced by a large number forcing the algorithm to backtrack,
i.e. decrease the linesearch factor or decrease the trust region radius.
When the algorithm is evaluating a numerical jacobian, it will stop with an error message
on detecting non-finite function values.

The \code{control} argument is a list that can supply any of the
following components:
\describe{
\item{\code{xtol}}{The relative steplength tolerance.
When the relative  steplength of all scaled x values is smaller than this value
convergence is declared. The default value is \eqn{10^{-8}}{1e-8}.
}
\item{\code{ftol}}{The function value tolerance.
Convergence is declared when the largest absolute function value is smaller than \code{ftol}.
The default value is \eqn{10^{-8}}{1e-8}.
}
\item{\code{btol}}{The backtracking tolerance.
When the relative steplength in a backtracking step to find an acceptable point is smaller
than the backtracking tolerance, the backtracking is terminated.
In the \code{Broyden} method a new Jacobian will be calculated if the Jacobian is outdated.
The default value is \eqn{10^{-3}}{1e-3}.
}
\item{\code{sigma}}{Reduction factor for the geometric linesearch. The default value is 0.5.}
\item{\code{scalex}}{a vector of scaling values for the parameters.
The inverse of a scale value is an indication of the size of a parameter.
The default value is 1.0 for all scale values.}
\item{\code{maxit}}{The maximum number of major iterations. The default value is 150.}
\item{\code{trace}}{Non-negative integer. A value of 1 will give a detailed report of the
progress of the iteration.}
\item{\code{chkjac}}{A logical value indicating whether to check an analytical jacobian, if
supplied. The default value is \code{FALSE}. The first 10 errors are printed.
The code for this check is derived from the code in Bouaricha and Schnabel (1997).}
\item{\code{delta}}{Initial (scaled) trust region radius.
A value of \eqn{-1.0} is replaced by the length of the Cauchy step in the initial point.
A value of \eqn{-2.0} is replaced by the length of the Newton step in the initial point.
Any value less than or equal to 0 and not equal to \eqn{-2.0}, will be replaced by \eqn{-1.0};
the algorithm will thus start with the length of the Cauchy step in the initial point.
If it is positive it will be set to the smaller of the value supplied or the maximum stepsize.
The default is \eqn{-2.0}.}
\item{\code{stepmax}}{Maximum scaled stepsize.
If this is negative then the maximum stepsize is set to the largest positive representable number.
The default is \eqn{-1.0}, so there is no default maximum stepsize.}
}
}
\value{
A list containing components
  \item{x}{final values for x}
  \item{fvec}{function values}
  \item{termcd}{termination code as integer.
         The values returned are
         \describe{
         \item{\code{1}}{Function criterion is near zero.
         Convergence of function values has been achieved.}
         \item{\code{2}}{x-values within tolerance. This means that the relative distance between two
         consecutive x-values is smaller than \code{xtol}.}
         \item{\code{3}}{No better point found.
         This means that the algorithm has stalled and cannot find an acceptable new point.
         This may or may not indicate acceptably small function values.}
         \item{\code{4}}{Iteration limit \code{maxit} exceeded.}
         \item{\code{5}}{Jacobian is too ill-conditioned.}
         \item{\code{6}}{Jacobian is singular.}
         }
   }
  \item{message}{a string describing the termination code}
  \item{scalex}{a vector containing the scaling factors,
                which will be the final values when automatic scaling was selected}
  \item{njcnt}{number of jacobian evaluations}
  \item{nfcnt}{number of function evaluations, excluding those required for calculating a jacobian}
}

%\note{
%}
%
\section{Warning}{You cannot use this function recursively.
Thus function \code{fn} should not in its turn call \code{nleqslv}.
}
\references{
Bouaricha, A. and Schnabel, R.B. (1997),
Algorithm 768: TENSOLVE: A Software Package for
                 Solving Systems of Nonlinear Equations and Nonlinear
                 Least-squares Problems Using Tensor Methods,
\emph{Transactions on Mathematical Software}, \bold{23}, 2, pp. 174--195.

Dennis, J.E. Jr  and Schnabel, R.B. (1996), \emph{Numerical Methods for Unconstrained Optimization
and Nonlinear Equations}, Siam.

\enc{Moré}{More}, J.J. (1978), The Levenberg-Marquardt Algorithm, Implementation and Theory,
In \emph{Numerical Analysis}, G.A. Watson (Ed.),
Lecture Notes in Mathematics 630, Springer-Verlag, pp. 105--116.

%%Nocedal, J. and Wright, S.J. (2006), \emph{Numerical Optimization}, Springer.

Powell, M.J.D. (1970), A hybrid method for nonlinear algebraic equations,
In \emph{Numerical Methods for Nonlinear Algebraic Equations}, P. Rabinowitz (Ed.), Gordon \& Breach.

Powell, M.J.D. (1970), A Fortran subroutine for solving systems nonlinear equations,
In \emph{Numerical Methods for Nonlinear Algebraic Equations}, P. Rabinowitz (Ed.), Gordon \& Breach.

Reichel, L. and W.B. Gragg (1990), Algorithm 686: FORTRAN subroutines for updating the QR decomposition,
\emph{ACM Trans. Math. Softw.}, \bold{16}, 4, pp. 369--377.
}
%
%\seealso{\code{\link{"BB"}}}
%
\examples{
# Dennis Schnabel example 6.5.1 page 149
dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}

jacdsln <- function(x) {
    n <- length(x)
    Df <- matrix(numeric(n*n),n,n)
    Df[1,1] <- 2*x[1]
    Df[1,2] <- 2*x[2]
    Df[2,1] <- exp(x[1]-1)
    Df[2,2] <- 3*x[2]^2

    Df
}

BADjacdsln <- function(x) {
    n <- length(x)
    Df <- matrix(numeric(n*n),n,n)
    Df[1,1] <- 4*x[1]
    Df[1,2] <- 2*x[2]
    Df[2,1] <- exp(x[1]-1)
    Df[2,2] <- 5*x[2]^2

    Df
}

xstart <- c(2,0.5)
fstart <- dslnex(xstart)
xstart
fstart

# a solution is c(1,1)

nleqslv(xstart, dslnex, control=list(btol=.01))

# Cauchy start
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta=-1.0))

# Newton start
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta=-2.0))

\dontrun{
# this gives a bus error on a ppc version of Mac OS X
# I have no opportunity to run ppc code so I can't debug

# example for error message in case of an error in analytical jacobian

nleqslv(xstart, dslnex, BADjacdsln, method="Newton",
                    control=list(btol=.01,chkjac=TRUE))
}

\dontrun{
# singular function with a singular start

sngsrt <- function(x) {
    y <- numeric(length(x))
    y[1] <- x[1]^2/2 + x[2]
    y[2] <- x[1] + x[2]^2/2

    y
}

sngsrtjac <- function(x)  {
    n <- length(x)
    Df <- matrix(numeric(n*n),n,n)
    Df[1,1] <- x[1]
    Df[1,2] <- 1
    Df[2,1] <- 1
    Df[2,2] <- x[2]

    Df
}
xstart <- c(1,1)
nleqslv(xstart, sngsrt, sngsrtjac, global="dbldog",
        control=list(trace=1,btol=.01))
}        
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonlinear}
\keyword{optimize}

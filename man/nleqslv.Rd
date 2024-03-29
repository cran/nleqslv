\name{nleqslv}
\alias{nleqslv}
\encoding{UTF-8}
\title{Solving systems of nonlinear equations with Broyden or Newton}
\description{
The function solves a system of nonlinear equations with either a Broyden or a
full Newton method. It provides line search and trust region global strategies
for difficult systems.
}
\usage{
nleqslv(x, fn, jac=NULL, ...,
               method = c("Broyden", "Newton"),
               global = c("dbldog", "pwldog",
                          "cline", "qline", "gline", "hook", "none"),
               xscalm = c("fixed","auto"),
               jacobian=FALSE,
               control = list()
	)
}
\arguments{
 \item{x}{A numeric vector with an initial guess of the root of the function.}
 \item{fn}{A function of \code{x} returning a vector of function values
           with the same length as the vector \code{x}.}
 \item{jac}{A function to return the Jacobian for the \code{fn} function.
            For a vector valued function \code{fn} the Jacobian must be a numeric
            matrix of the correct dimensions.
            For a scalar valued function \code{fn} the \code{jac} function may return a scalar.
            If not supplied numerical derivatives will be used.}
 \item{\dots}{Further arguments to be passed to \code{fn} and \code{jac}.}
 \item{method}{The method to use for finding a solution. See \sQuote{Details}.}
 \item{global}{The global strategy to apply. See \sQuote{Details}.}
 \item{xscalm}{The type of x scaling to use. See \sQuote{Details}.}
 \item{jacobian}{A logical indicating if the estimated (approximate) jacobian in the solution should be returned.
 See \sQuote{Details}.}
 \item{control}{A named list of control parameters. See \sQuote{Details}.}
}
\details{
The algorithms implemented in \code{nleqslv} are based on Dennis and Schnabel (1996).

\subsection{Methods}{

Method \code{Broyden} starts with a computed Jacobian of the function and
then updates this Jacobian after each successful iteration using the so-called Broyden update.
This method often shows super linear convergence towards a solution.
When \code{nleqslv} determines that it cannot
continue with the current Broyden matrix it will compute a new Jacobian.

Method \code{Newton} calculates a Jacobian of the function \code{fn} at each iteration.
Close to a solution this method will usually show quadratic convergence.

Both methods apply a so-called (backtracking) global strategy to find
a better (more acceptable) iterate.  The function criterion used by the algorithm is half of
the sum of squares of the function values and \dQuote{acceptable} means sufficient decrease
of the current function criterion value compared to that of the previous iteration.
A comprehensive discussion of these issues can be found in Dennis and Schnabel (1996).
Both methods apply an unpivoted QR-decomposition to the Jacobian as implemented in Lapack.
The Broyden method applies a rank-1 update to the Jacobian at the end of each iteration
and is based on a simplified and modernized version of the algorithm described in Reichel and Gragg (1990).
}

\subsection{Global strategies}{

When applying a full Newton or Broyden step does not yield a sufficiently smaller
function criterion value \code{nleqslv} will attempt to decrease the steplength using
one of several so-called global strategies.

The \code{global} argument indicates which global strategy to use or to use no global strategy
\describe{
    \item{\code{cline}}{a cubic line search}
    \item{\code{qline}}{a quadratic line search}
    \item{\code{gline}}{a geometric line search}
    \item{\code{dbldog}}{a trust region method using the double dogleg method as described in
                         Dennis and Schnabel (1996)}
    \item{\code{pwldog}}{a trust region method using the Powell dogleg method
                         as developed by Powell (1970).}
    \item{\code{hook}}{a trust region method described by Dennis and Schnabel (1996) as
                       \emph{The locally constrained optimal (\dQuote{hook}) step}.
                       It is equivalent to a Levenberg-Marquardt algorithm as described in
                       \enc{Moré}{More} (1978) and Nocedal and Wright (2006).}
    \item{\code{none}}{Only a pure local Newton or Broyden iteration is used.
                       The maximum stepsize (see below) is taken into account.
                       The default maximum number of iterations (see below) is set to 20.}
}
The double dogleg method is the default global strategy employed by this package.

Which global strategy to use in a particular situation is a matter of trial and error.
When one of the trust region methods fails, one of the line search strategies should be tried.
Sometimes a trust region will work and sometimes a line search method; neither has a clear advantage but
in many cases the double dogleg method works quite well.

When the function to be solved returns non-finite function values for a parameter vector \code{x}
and the algorithm is \emph{not} evaluating a numerical Jacobian, then any non-finite
values will be replaced by a large number forcing the algorithm to backtrack,
i.e. decrease the line search factor or decrease the trust region radius.
}

\subsection{Scaling}{

The elements of vector \code{x} may be scaled during the search for a zero of \code{fn}.
The  \code{xscalm} argument provides two possibilities for scaling
\describe{
    \item{\code{fixed}}{the scaling factors are set to the values supplied in
    the \code{control} argument and remain unchanged during the iterations.
    The scaling factor of any element of \code{x}  should be set
    to the inverse of the typical value of that element of \code{x},
    ensuring that all elements of \code{x} are approximately equal in size.}

    \item{\code{auto}}{the scaling factors
    are calculated from the euclidean norms of the columns of the Jacobian matrix.
    When a new Jacobian is computed, the scaling values will be set to the
    euclidean norm of the corresponding column if that is larger than the current scaling value.
    Thus the scaling values will not decrease during the iteration.
    This is the method described in \enc{Moré}{More} (1978).
    Usually manual scaling is preferable.}
}
}

\subsection{Jacobian}{

When evaluating a numerical Jacobian, an error message will be issued
on detecting non-finite function values.
An error message will also be issued when a user supplied jacobian contains non-finite entries.

When the \code{jacobian} argument is set to \code{TRUE} the final Jacobian or Broyden matrix
will be returned in the return list.
The default value is \code{FALSE}; i.e. to not return the final matrix.
There is no guarantee that the final Broyden matrix resembles the actual Jacobian.

The package can cope with a singular or ill-conditioned Jacobian if needed
by setting the \code{allowSingular} component of the \code{control} argument.
The method used is described in Dennis and Schnabel (1996);
it is equivalent to a Levenberg-Marquardt type adjustment with a small damping factor.
\emph{There is no guarantee that this method will be successful.}
Warning: \emph{\code{nleqslv} may report spurious convergence in this case.}

By default \code{nleqslv} returns an error
if a Jacobian becomes singular or very ill-conditioned.
A Jacobian is considered to be very ill-conditioned when the estimated inverse
condition is less than or equal to a specified tolerance with a default value
equal to \eqn{10^{-12}}{1e-12}; this can be changed and made smaller
with the \code{cndtol} item of the \code{control} argument.
\emph{There is no guarantee that any change will be effective.}
}

\subsection{Control options}{

The \code{control} argument is a named list that can supply any of the
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

    \item{\code{cndtol}}{The tolerance of the test for ill conditioning of
    the Jacobian or Broyden approximation. If less than the machine precision it will
    be silently set to the machine precision.
    When the estimated inverse condition of the (approximated) Jacobian matrix is less than or equal to
    the value of \code{cndtol} the matrix is deemed to be ill-conditioned,
    in which case an error will be reported if the \code{allowSingular} component is set to \code{FALSE}.
    The default value is \eqn{10^{-12}}{1e-12}.
    }

    \item{\code{sigma}}{Reduction factor for the geometric line search. The default value is 0.5.}

    \item{\code{scalex}}{a vector of scaling values for the parameters.
    The inverse of a scale value is an indication of the size of a parameter.
    The default value is 1.0 for all scale values.}

    \item{\code{maxit}}{The maximum number of major iterations.
    The default value is 150 if a global strategy has been specified.
    If no global strategy has been specified the default is 20.}

    \item{\code{trace}}{Non-negative integer. A value of 1 will give a detailed report of the
    progress of the iteration. For a description see \code{\link{Iteration report}}.}

    \item{\code{chkjac}}{A logical value indicating whether to check a user supplied Jacobian, if
    supplied. The default value is \code{FALSE}. The first 10 errors are printed.
    The code for this check is derived from the code in Bouaricha and Schnabel (1997).}

    \item{\code{delta}}{Initial (scaled) trust region radius.
    A value of \eqn{-1.0} or \code{"cauchy"} is replaced by the length of the Cauchy step in the initial point.
    A value of \eqn{-2.0} or \code{"newton"} is replaced by the length of the Newton step in the initial point.
    Any numeric value less than or equal to 0 and not equal to \eqn{-2.0}, will be replaced by \eqn{-1.0};
    the algorithm will then start with the length of the Cauchy step in the initial point.
    If it is numeric and positive it will be set to the smaller of the value supplied or the maximum stepsize.
    If it is not numeric and not one of the permitted character strings then an error message will be issued.
    The default is \eqn{-2.0}.}

    \item{\code{stepmax}}{Maximum scaled stepsize.
    If this is negative then the maximum stepsize is set to the largest positive representable number.
    The default is \eqn{-1.0}, so there is no default maximum stepsize.}

    \item{\code{dsub}}{Number of non zero subdiagonals of a banded Jacobian.
    The default is to assume that the Jacobian is \emph{not} banded.
    Must be specified if \code{dsuper} has been specified and must be larger than zero when \code{dsuper} is zero.}

    \item{\code{dsuper}}{Number of non zero super diagonals of a banded Jacobian.
    The default is to assume that the Jacobian is \emph{not} banded.
    Must be specified if \code{dsub} has been specified and must be larger than zero when \code{dsub} is zero.}

    \item{\code{allowSingular}}{A logical value indicating if a small correction
    to the Jacobian when it is singular or too ill-conditioned is allowed.
    If the correction is less than \code{100*.Machine$double.eps} the correction
    cannot be applied and an unusable Jacobian will be reported.
    The method used is similar to a Levenberg-Marquardt correction and
    is explained in Dennis and Schnabel (1996) on page 151.
    It may be necessary to choose a higher value for \code{cndtol} to enforce the correction.
    The default value is \code{FALSE}.
    }
}
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
         consecutive x-values is smaller than \code{xtol} but that
         the function value criterion is still larger than \code{ftol}.
         \emph{Function values may not be near zero; therefore the user must check if
               function values are acceptably small.}}
         \item{\code{3}}{No better point found.
         This means that the algorithm has stalled and cannot find an acceptable new point.
         This may or may not indicate acceptably small function values.}
         \item{\code{4}}{Iteration limit \code{maxit} exceeded.}
         \item{\code{5}}{Jacobian is too ill-conditioned.}
         \item{\code{6}}{Jacobian is singular.}
         \item{\code{7}}{Jacobian is unusable.}
         \item{\code{-10}}{User supplied Jacobian is most likely incorrect.}
         }
    }
    \item{message}{a string describing the termination code}
    \item{scalex}{a vector containing the scaling factors,
                which will be the final values when automatic scaling was selected}
    \item{njcnt}{number of Jacobian evaluations}
    \item{nfcnt}{number of function evaluations, excluding those required for calculating a Jacobian and
               excluding the initial function evaluation (at iteration 0)}
    \item{iter}{number of outer iterations used by the algorithm. This excludes the initial iteration.
    The number of backtracks can be calculated as the difference between the \code{nfcnt} and \code{iter}
    components.}
    \item{jac}{the final Jacobian or the Broyden approximation if \code{jacobian} was set to \code{TRUE}.
    If no iterations were executed, as may happen when the initial guess are sufficiently close
    the solution, there is no Broyden approximation and the returned matrix will always be the actual Jacobian.
    If the matrix is singular or too ill-conditioned the returned matrix is of no value.}
}

%\note{
%}
%
\section{Warning}{You cannot use this function recursively.
Thus function \code{fn} should not in its turn call \code{nleqslv}.
}

\seealso{If this function cannot solve the supplied function then it is
a good idea to try the function \link{testnslv} in this package.
For detecting multiple solutions see \link{searchZeros}.}

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

Golub, G.H and C.F. Van Loan (1996), Matrix Computations (3rd edition), The John Hopkins University Press.

Higham, N.J. (2002), Accuracy and Stability of Numerical Algorithms, 2nd ed., SIAM, pp. 10--11.

Nocedal, J. and Wright, S.J. (2006), \emph{Numerical Optimization}, Springer.

Powell, M.J.D. (1970), A hybrid method for nonlinear algebraic equations,
In \emph{Numerical Methods for Nonlinear Algebraic Equations}, P. Rabinowitz (Ed.), Gordon & Breach.

Powell, M.J.D. (1970), A Fortran subroutine for solving systems nonlinear equations,
In \emph{Numerical Methods for Nonlinear Algebraic Equations}, P. Rabinowitz (Ed.), Gordon & Breach.

% QRupdate: a Fortran library for fast updating of QR and Cholesky decompositions,
% \url{http://sourceforge.net/projects/qrupdate/}.

Reichel, L. and W.B. Gragg (1990), Algorithm 686: FORTRAN subroutines for updating the QR decomposition,
\emph{ACM Trans. Math. Softw.}, \bold{16}, 4, pp. 369--377.
}

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
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="cauchy"))

# Newton start
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))

# final Broyden approximation of Jacobian (quite good)
z <- nleqslv(xstart, dslnex, jacobian=TRUE,control=list(btol=.01))
z$x
z$jac
jacdsln(z$x)

# different initial start; not a very good final approximation
xstart <- c(0.5,2)
z <- nleqslv(xstart, dslnex, jacobian=TRUE,control=list(btol=.01))
z$x
z$jac
jacdsln(z$x)

\dontrun{
# no global strategy but limit stepsize
# but look carefully: a different solution is found
nleqslv(xstart, dslnex, method="Newton", global="none", control=list(trace=1,stepmax=5))

# but if the stepsize is limited even more the c(1,1) solution is found
nleqslv(xstart, dslnex, method="Newton", global="none", control=list(trace=1,stepmax=2))

# Broyden also finds the c(1,1) solution when the stepsize is limited
nleqslv(xstart, dslnex, jacdsln, method="Broyden", global="none", control=list(trace=1,stepmax=2))
}

# example with a singular jacobian in the initial guess
f <- function(x) {
    y <- numeric(3)
    y[1] <- x[1] + x[2] - x[1]*x[2] - 2
    y[2] <- x[1] + x[3] - x[1]*x[3] - 3
    y[3] <- x[2] + x[3] - 4
    return(y)
}

Jac <- function(x) {
    J <- matrix(0,nrow=3,ncol=3)
    J[,1] <- c(1-x[2],1-x[3],0)
    J[,2] <- c(1-x[1],0,1)
    J[,3] <- c(0,1-x[1],1)
    J
}

# exact solution
xsol <- c(-.5, 5/3 , 7/3)
xsol

xstart <- c(1,2,3)
J <- Jac(xstart)
J
rcond(J)

z <- nleqslv(xstart,f,Jac, method="Newton",control=list(trace=1,allowSingular=TRUE))
all.equal(z$x,xsol)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonlinear}
\keyword{optimize}

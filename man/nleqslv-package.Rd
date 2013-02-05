\name{nleqlsv-package}
\alias{nleqlsv-package}
\alias{nleqlsv.Intro}
\docType{package}

\title{Solving Nonlinear Systems of Equations}

\description{
The \pkg{nleqlsv} package provides two algorithms for solving (dense) non linear systems of equations.
The methods provided are
  \itemize{
    \item a Broyden Secant method where the matrix of derivatives is updated after each major iteration
          using the Broyden rank 1 update.
    \item a full Newton method where the Jacobian matrix of derivatives is recalculated at each iteration
  }
Both methods utilise global strategies such as line search or trust region methods
whenever the standard Newton/Broyden step does not lead to a point closer to a root
of the equation system.
Line search may be either quadratic or geometric. The trust region methods are either
the double dogleg or the Powell single dogleg method.

The algorithms provided in this package are derived from Dennis and Schnabel (1996).
The code is written in Fortran 77 and Fortran 95
and uses Lapack and BLAS routines as provided by the R system.
}

\author{Berend Hasselman \email{bhh@xs4all.nl}}

\references{
Dennis, J.E. Jr  and Schnabel, R.B. (1996), \emph{Numerical Methods for Unconstrained Optimization
and Nonlinear Equations}, Siam.
}
\keyword{package}

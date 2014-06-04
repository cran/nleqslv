
\name{print.test.nleqslv}
\title{Printing the result of \code{testnslv}}
\alias{print}
\alias{print.test.nleqslv}
\description{
  Print method for \code{test.nleqslv} objects.
}
\usage{
\method{print}{test.nleqslv}(x, digits=4, \dots)
}
\arguments{
  \item{x}{a \code{test.nleqslv} object}
  \item{digits}{specifies the minimum number of significant digits to be printed in values.}
  \item{\dots}{additional arguments to \code{print}.}
}
\details{
This is the \code{print} method for objects inheriting from
class \code{test.nleqslv}. It prints the call to \code{testnslv} followed by the description of the experiment
(if the \code{title} argument was specified in the call to \code{testnslv})
and the dataframe containing the results of \code{testnslv}.
}
\value{
It returns the object \code{x} invisibly.
}
\examples{
dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}
xstart <- c(1.5,0.5)
fstart <- dslnex(xstart)
z <- testnslv(xstart,dslnex)
print(z)
}
\keyword{print}

\name{nleqslv-iterationreport}
\alias{Iteration report}
\title{Detailed iteration report of nleqslv}
\description{ The format of the iteration report
provided by nleqslv when the \code{trace} component of the control parameter
has been set to 1.
}
\details{ 
All iteration reports consist of a series of columns with a header summarising the contents.
Common column headers are
\describe{
\item{\code{Iter}}{Iteration counter}
\item{\code{Jac}}{Jacobian type. The Jacobian type is indicated by \code{N} for a Newton Jacobian
or \code{B} for a Broyden updated matrix; optionally followed by the letter \code{s} 
indicating a totally singular matrix or the letter \code{i} indicating an ill-conditioned matrix.
Unless the Jacobian is singular, the Jacobian type is followed by an estimate
of the inverse condition number of the Jacobian in parentheses as computed by Lapack.
This column will be blank when backtracking is active.}
\item{\code{Fnorm}}{square of the euclidean norm of function values / 2}
\item{\code{Largest |f|}}{infinity norm of \eqn{f(x)}{f(x)} at the current point}
}
A sample iteration report for \code{qline} is
(some intercolumn space has been removed to make the table fit)  
\preformatted{
Iter        Jac Lambda        Ftarg         Fnorm   Largest |f|
   0                                 2.886812e+00  2.250000e+00
   1 N(9.6e-03) 1.0000 2.886235e+00  5.787362e+05  1.070841e+03
   1            0.1000 2.886754e+00  9.857947e+00  3.214799e+00
   1            0.0100 2.886806e+00  2.866321e+00  2.237878e+00
   2 B(2.2e-02) 1.0000 2.865748e+00  4.541965e+03  9.341610e+01
}
The column headed by \code{Lambda} shows the value of the line search parameter.
The column headed by \code{Ftarg} follows from a sufficient decrease requirement and
is the value below which \code{Fnorm} must drop if the current step is to be accepted.

The iteration report for \code{none} is almost the same as the above report, except that the column
labelled  \code{Ftarg} is omitted. The column \code{Lambda} gives the ratio of
the maximum stepsize and the length of the step taken by the algorithm.
It is either 1.0 exactly or smaller.

A sample iteration report for \code{dbldog} is 
(some intercolumn space has been removed to make the table fit)
\preformatted{
Iter        Jac   Lambda    Eta   Dlt0   Dltn        Fnorm  Largest |f|
   0                                          2.886812e+00 2.250000e+00
   1 N(9.6e-03) C        0.9544 0.4671 0.9343 1.699715e-01 5.421673e-01
   1            W 0.0833 0.9544 0.9343 0.4671 1.699715e-01 5.421673e-01
   2 B(1.1e-02) W 0.1154 0.4851 0.4671 0.4671 1.277667e-01 5.043571e-01
   3 B(7.3e-02) W 0.7879 0.7289 0.4671 0.0759 5.067893e-01 7.973542e-01
   3            C        0.7289 0.0759 0.1519 5.440250e-02 2.726084e-01
   4 B(8.3e-02) W 0.5307 0.3271 0.1519 0.3037 3.576547e-02 2.657553e-01
   5 B(1.8e-01) N        0.6674 0.2191 0.4383 6.566182e-03 8.555110e-02
   6 B(1.8e-01) N        0.9801 0.0376 0.0752 4.921645e-04 3.094104e-02
   7 B(1.9e-01) N        0.7981 0.0157 0.0313 4.960629e-06 2.826064e-03
   8 B(1.6e-01) N        0.3942 0.0029 0.0058 1.545503e-08 1.757498e-04
   9 B(1.5e-01) N        0.6536 0.0001 0.0003 2.968676e-11 5.983765e-06
  10 B(1.5e-01) N        0.4730 0.0000 0.0000 4.741792e-14 2.198380e-07
  11 B(1.5e-01) N        0.9787 0.0000 0.0000 6.451792e-19 8.118586e-10
}

After the column for the Jacobian the letters indicate the following
\describe{
\item{\code{C}}{a fraction (\eqn{\le1.0}{<=1.0}) of the Cauchy or steepest descent step is taken where
the fraction is the ratio of the trust region radius and the Cauchy steplength.}
\item{\code{W}}{a convex combination of the Cauchy and \code{eta}*(Newton step) is taken.
The number in the column headed by \code{Lambda} is the weight of the partial
Newton step.}
\item{\code{P}}{a fraction (\eqn{\ge1.0}{>=1.0}) of the partial Newton step, equal to \code{eta}*(Newton step),
is taken where
the fraction is the ratio of the trust region radius and the partial Newton steplength.} 
\item{\code{N}}{a normal full Newton step is taken.}
}
The number in the  column headed by \code{Eta} is calculated from an upper limit on the ratio
of the length of the steepest descent direction and the length of the Newton step.
See Dennis and Schnabel (1996) pp.139ff for the details.

The column headed by \code{Dlt0} gives the size of the trust region at the start of the current
iteration. 
The column headed by \code{Dltn} gives the size of the trust region when the current
step has been accepted by the dogleg algorithm. 
The size of the trust region is decreased when the actual reduction of the function value norm
does not agree well with the predicted reduction from the linear approximation of the function.
And increased when the actual and predicted reduction are sufficiently close.

Normally the initial trust region radius is the same as 
the final trust region radius of the previous iteration but the size of the trust region is
restricted by the size of the current Newton step. So when full Newton steps
are being taken, the trust region radius at the start of an iteration may be less than the
final value of the previous iteration.

A sample iteration report for \code{pwldog} is
(some intercolumn space has been removed to make the table fit)  
\preformatted{
Iter        Jac   Lambda   Dlt0   Dltn        Fnorm  Largest |f|
   0                                   2.886812e+00 2.250000e+00
   1 N(9.6e-03) C        0.4671 0.9343 1.699715e-01 5.421673e-01
   1            W 0.0794 0.9343 0.4671 1.699715e-01 5.421673e-01
   2 B(1.1e-02) W 0.0559 0.4671 0.4671 1.205661e-01 4.890487e-01
   3 B(7.3e-02) W 0.5662 0.4671 0.0960 4.119560e-01 7.254441e-01
   3            W 0.0237 0.0960 0.1921 4.426507e-02 2.139252e-01
   4 B(8.8e-02) W 0.2306 0.1921 0.3842 2.303135e-02 2.143943e-01
   4            W 0.4769 0.3842 0.1921 2.303135e-02 2.143943e-01
   5 B(1.9e-01) N        0.1375 0.2750 8.014508e-04 3.681498e-02
}
This is much simpler than the double dogleg report, since the (single) dogleg takes
either a steepest descent step, a convex combination of the steepest descent and Newton steps
or a full Newton step. The number in the column \code{Lambda} is the weight of the Newton step.
It is a special case of the double dogleg method with \code{eta} equal to 1.
}
\seealso{\code{\link{nleqslv}}}

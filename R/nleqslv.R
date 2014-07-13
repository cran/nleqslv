#
# Interface to Fortran library for solving a system of non linear equations
# with either a Broyden or a full Newton method
# There a four global search methods:
#   quadratic linesearch, geometric linesearch
#   double dogleg trust region a la Dennis Schnabel
#   powell single dogleg a la Minpack
#
#

nleqslv <- function(x, fn, jac = NULL, ...,
                    method = c("Broyden", "Newton"),
                    global = c("dbldog", "pwldog", "cline", "qline", "gline", "none"),
                    xscalm = c("fixed","auto"),
                    jacobian=FALSE,
                    control = list())
{
    fn1  <- function(par) fn(par, ...)
    jac1 <- if (!is.null(jac)) function(par) jac(par, ...)

    method <- match.arg(method)
    global <- match.arg(global)
    xscalm <- match.arg(xscalm)

    ## Defaults
    con <- list(ftol=1e-8, xtol=1e-8,
                btol=1e-3,
                stepmax=-1.0, delta=-2.0, sigma=0.5,
                scalex = rep(1, length(x)),
                maxit=150,
                trace=0,
                chkjac=FALSE,
                cndtol=1e-12,
                dsub=-1L,
                dsuper=-1L
               )

    # limit maximum number of iterations for pure local strategy
    if( global == "none" ) con$maxit=20

    # check names of control argument
    namc <- names(control)
    if (!all(namc %in% names(con)))
        stop("unknown names in control: ", paste(namc[!(namc %in% names(con))], collapse=", "))
    con[namc] <- control

    tmp <- con[["delta"]]
    if( is.character(tmp) ) {
        k <- match(tolower(tmp), c("cauchy","newton"))
        con[["delta"]] <- as.numeric(-k)
    }
    else if( !is.numeric(tmp) ) stop('control["delta"] should be either character or numeric')

    # to reset flag for checking recursive calls (not allowed for now)
    on.exit(.C("deactivatenleq",PACKAGE="nleqslv"))
    out <- .Call("nleqslv", x, fn1, jac1, method, global, xscalm, jacobian, con, new.env(), PACKAGE = "nleqslv")

    out
}

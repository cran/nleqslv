#
# Interface to Fortran library for solving a system of non linear equations
# with either a Broyden or a full Newton method
# There a four global search methods:
#   quadratic linesearch, geometric linesearch
#   double dogleg trust region a la Dennis Schnabel
#   powell single dogleg a la Minpack
#
#

.onLoad <- function(lib,pkg) {
    library.dynam("nleqslv","nleqslv")
}

nleqslv <- function(x, fn, jac = NULL, ...,
                    method = c("Broyden", "Newton"),
                    global = c("dbldog", "pwldog", "qline", "gline"),
                    xscalm = c("fixed","auto"),
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
                chkjac=FALSE
               )                     
               
    # check names of control argument
    namc <- names(control)
    if (!all(namc %in% names(con))) 
        stop("unknown names in control: ", paste(namc[!(namc %in% names(con))], collapse=", "))
    con[namc] <- control

    # to reset flag for checking recursive calls (not allowed for now)
    on.exit(.C("deactivatenleq",PACKAGE="nleqslv"))
    out <- .Call("nleqslv", x, fn1, jac1, method, global, xscalm, con, new.env(), PACKAGE = "nleqslv")

    out
}

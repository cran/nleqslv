
testnslv <- function(x, fn, jac=NULL, ...,
                            method=c("Newton", "Broyden"),
                            global=c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"),
                            Nrep=0L, title=NULL
                          )
{
    # utility functions
    catmsg <- function(m,g,res) {
        cat(sprintf("Error (method=%s global=%s): %s\n",m,g,attr(res,"condition")$message))
    }

    makeerrlist <- function(m,g,cpusecs=NULL) {
        if(is.null(cpusecs)) {
            z <- list(method=m, global=g, termcd=NA, fctcall=NA, jaccall=NA, itercnt=NA, message="ERROR",fnorm=NA)
        } else {
            z <- list(method=m, global=g, termcd=NA, fctcall=NA, jaccall=NA, itercnt=NA, message="ERROR",fnorm=NA,
                             cpusecs=cpusecs)
        }
        z
    }

    makereslist <- function(m,g,res,cpusecs=NULL) {
        fnorm <- sum(res$fvec^2)/2
        if(res$termcd < 0 ) stop("User supplied jacobian most likely incorrect: cannot continue") else
            message <- switch(res$termcd, "Fcrit", "Xcrit", "Stalled", "Maxiter", "Illcond", "Singular")
        
        if(is.null(cpusecs)) {
           z <- list(method=m, global=g, termcd=res$termcd, fctcall=res$nfcnt, jaccall=res$njcnt, itercnt=res$iter,
                            message=message, fnorm=fnorm)

        } else {
           z <- list(method=m, global=g, termcd=res$termcd, fctcall=res$nfcnt, jaccall=res$njcnt, itercnt=res$iter,
                            message=message, fnorm=fnorm, cpusecs=cpusecs)
        }

    }

    methods <- match.arg(method, c("Newton", "Broyden"), several.ok=TRUE)
    globals <- match.arg(global, c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"), several.ok=TRUE)

    my.call <- match.call()
    reslist <- vector("list", length(methods)*length(globals))

    idx <- 1
    for(m in methods)
        for(g in globals) {
            if( Nrep >= 1) {
                mytime <- system.time( for(k in seq_len(Nrep)) {
                                res <- try(nleqslv(x, fn, jac, ..., method=m, global=g), silent=TRUE)
                                if(inherits(res,"try-error")) break
                            }, gcFirst = FALSE)
                cpus <- mytime[3]
            } else {
                res <- try(nleqslv(x, fn, jac, ..., method=m, global=g),silent=TRUE)
                cpus <- NULL
            }
            if(inherits(res,"try-error")) {
                catmsg(m,g,res)
                z <- makeerrlist(m,g,cpus)
            } else {
                z <- makereslist(m,g,res,cpus)
            }
            reslist[[idx]] <- z
            idx <- idx+1
        }

# from http://stackoverflow.com/questions/4512465/what-is-the-most-efficient-way-to-cast-a-list-as-a-data-frame?rq=1

    ## @Martin Morgan's Map() sapply() solution:
    f <- function(x) function(i) sapply(x, `[[`, i)
    z <- as.data.frame(Map(f(reslist), names(reslist[[1]])), stringsAsFactors=FALSE)

    res <- list()
    res$out <- z
    res$call <- my.call
    res$title <- title
    class(res) <- "test.nleqslv"
    res
}


print.test.nleqslv <- function(x, digits=4, ...) {
    if(!inherits(x, "test.nleqslv"))
        stop("method is only for test.nleqslv objects")

    cat("Call:\n",paste0(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
    if(is.null(x$title)) cat("Results:\n") else cat("Results: ",x$title,"\n", sep="")
    print(x$out, digits=digits,...)
    invisible(x)
}

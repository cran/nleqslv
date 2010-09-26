# Dennis Schnabel example

library("nleqslv")
    
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

xstart <- c(2,0.5)
fstart <- dslnex(xstart)
xstart
fstart

# a solution is c(1,1)
znleqa <- nleqslv(xstart, dslnex, global="qline",xscalm="auto",control=list(trace=1,btol=.01))
znleqa

znleqb <- nleqslv(xstart, dslnex, global="gline",xscalm="auto",control=list(trace=1,btol=.01))
znleqb

znleqc <- nleqslv(xstart, dslnex, global="dbldog",xscalm="auto",
                    control=list(trace=1,btol=.01,delta=-1.0))
znleqc

znleqd <- nleqslv(xstart, dslnex, global="pwldog",xscalm="auto",
                    control=list(trace=1,btol=.01,delta=-1.0))
znleqd

znleqe <- nleqslv(xstart, dslnex, global="dbldog",xscalm="auto",
                    control=list(trace=1,btol=.01,delta=-2.0))
znleqe

znleqf <- nleqslv(xstart, dslnex, global="pwldog",xscalm="auto",
                    control=list(trace=1,btol=.01,delta=-2.0))
znleqf

znlejc <- nleqslv(xstart, dslnex, jacdsln, global="dbldog",xscalm="auto",
                    control=list(trace=1,btol=.01,delta=-1.0))
znlejc

znlejd <- nleqslv(xstart, dslnex, jacdsln, global="dbldog",xscalm="auto",
                    method="Newton", control=list(trace=1,btol=.01,delta=-1.0))
znlejd

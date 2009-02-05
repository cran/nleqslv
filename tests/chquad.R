
# Chebyquad functions (no solution for n=8)

library("nleqslv")

chebyquad <- function(x) {
	n <- length(x)
	y <- numeric(n)
 
  	for(j in 1:n) {
    	t1 <- 1.0
    	t2 <- 2.0*x[j] - 1.0
        tmp <- 2.0*t2
        
    	for(i in 1:n) {
      		y[i] <- y[i] + t2
      		t3 <- tmp * t2 - t1
      		t1 <- t2
      		t2 <- t3
    	}
  	}

  	y = y / n

  	for(i in 1:n) {
    	if ( i%%2 == 0 ) {
      		y[i] = y[i] + 1.0 / (i * i - 1)
    	}
  	}  

	y
}

chebyinit <- function(n) {
    x <- (1:n) / (n + 1)
}

for (k in 1:9) {
	n <- k
	if( k != 8 ) {
		xstart <- chebyinit(n)
		fstart <- chebyquad(xstart)

		zz <- nleqslv(xstart, chebyquad, global="dbldog",
		               control=list(ftol=1e-8,xtol=1e-8,trace=1,btol=.01,delta=-2))
		print(zz)
	}
}

#include <R.h>
#include <Rdefines.h>

typedef struct opt_struct {
    SEXP x;
    SEXP fcall;
    SEXP jcall;
    SEXP env;
} opt_struct, *OptStruct;

OptStruct OS;

/*
 * set to 1 if main function called
 * used to check for a recursive call
 * which is not allowed (for the time being)
 */
static int activeflag = 0;

void deactivatenleq(void);  /* called on.exit; see .R code */

void fcnval(double *xc, double *fc, int *n, int *flag);
void fcnjac(double *rjac, int *ldr, double *x, int *n);

void F77_NAME(nwnleq)(double *x, int *n, double *scalex, int *maxit, int *jacflg, double *xtol,
                      double *ftol, double *btol,
                      int *method, int *global, int *xscalm, double *stepmx, double *dlt, double *sigma,
                      double *rwork, int *lrwork, double *rcdwrk, int *icdwrk, double *qrwork, int *qrwsiz,
                      void (*fcnjac)(double *rjac, int *ldr, double *x, int *n),
                      void (*fcnval)(double *xc, double *fc, int *n, int *flag),
   			          int *outopt, double *xp, double *fp, double *gp, int *njcnt, int *nfcnt, int *termcd);

void F77_NAME(nwqmem)(int *n, int *wrksiz); /* returns size of optimal QR work memory */

void deactivatenleq(void)
{
    activeflag = 0;
}

static SEXP getListElement(SEXP list, char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    	elmt = VECTOR_ELT(list, i);
	    	break;
		}
    return elmt;
}

static double *real_vector(int n)
{
    return (double*) R_alloc(n, sizeof(double));
}

static int *int_vector(int n)
{
    return (int*) R_alloc(n, sizeof(int));
}

static void trace_header(int method, int global, int xscalm, double sigma,
                         double dlt, double stepmx, double ftol, double xtol, double btol)
{	/* print header for iteration tracing output for documentation purposes
	 */

	Rprintf("  Algorithm parameters\n  --------------------\n");
	Rprintf("  Method: %s", method == 1 ? "Broyden" : "Newton");
	Rprintf("  Global strategy: ");
	switch(global)
	{
        case 0: Rprintf("quadratic linesearch\n"); break;
	    case 1: Rprintf("geometric linesearch (reduction = %g)\n", sigma); break;
	    case 2: Rprintf("double dogleg (initial trust region = %g)\n", dlt); break;
	    case 3: Rprintf("single dogleg (initial trust region = %g)\n", dlt); break;
        default: error("Internal: invalid global value in trace_header\n");
    }

	Rprintf("  Maximum stepsize = %g\n", stepmx);
    Rprintf("  Scaling: %s\n", xscalm == 0 ? "fixed" : "automatic");

	Rprintf("  ftol = %g xtol = %g btol = %g\n\n", ftol,xtol,btol);
    Rprintf("  Iteration report\n  ----------------\n");
}

static char *fcn_message(char *msg, int termcd)
{
    if      (termcd == 1)
        sprintf(msg, "Function criterion near zero");
    else if (termcd == 2)
        sprintf(msg, "x-values within tolerance `xtol'");
    else if (termcd == 3)
        sprintf(msg, "No better point found (algorithm has stalled)");
    else if (termcd == 4)
        sprintf(msg, "Iteration limit exceeded");
    else if (termcd == -10 )
        sprintf(msg, "Analytical Jacobian most likely incorrect");
    else
        sprintf(msg, "`termcd' == %d should *NEVER* be returned! Please report bug to <bhh@xs4all.nl>.", termcd);

    return msg;
}

/*
 * interface to user supplied R function
 */

void fcnval(double *xc, double *fc, int *n, int *flag)
{
    int i;
    SEXP sexp_fvec;

    if (IS_NUMERIC(OS->x))
        for (i = 0; i < *n; i++) {
            NUMERIC_POINTER(OS->x)[i] = xc[i];
    }
    else
        for (i = 0; i < *n; i++) {
            NUMERIC_POINTER(VECTOR_ELT(OS->x, i))[0] = xc[i];
    }

    SETCADR(OS->fcall, OS->x);
    PROTECT(sexp_fvec = eval(OS->fcall, OS->env));

    for (i = 0; i < *n; i++) {
        fc[i] = NUMERIC_POINTER(sexp_fvec)[i];
		if( !R_FINITE(fc[i]) ) {
            fc[i] = sqrt(DBL_MAX / (double)(*n)); /* should force backtracking */
		    if( *flag ) {
				error("Non-finite value(s) detected in jacobian (row=%d,col=%d)",i+1,*flag);
	        }
        }
	}

    UNPROTECT(1);
}

void FCNJACDUM(double *rjac, int *ldr, double *x, int *n)
{
	error("FCNJACDUM should not have been called");
}

/*
 * interface to user supplied jacobian function
*/

void fcnjac(double *rjac, int *ldr, double *x, int *n)
{
    int i, j;
    SEXP sexp_fjac;

    if (IS_NUMERIC(OS->x))
        for (i = 0; i < *n; i++) {
            if (!R_FINITE(x[i]))
                error("non-finite value supplied by Nwnleq!");
            NUMERIC_POINTER(OS->x)[i] = x[i];
        }
    else
        for (i = 0; i < *n; i++) {
            if (!R_FINITE(x[i]))
                error("non-finite value supplied by Nwnleq!");
            NUMERIC_POINTER(VECTOR_ELT(OS->x, i))[0] = x[i];
        }

    SETCADR(OS->jcall, OS->x);
    PROTECT(sexp_fjac = eval(OS->jcall, OS->env));

    for (j = 0; j < *n; j++)
        for (i = 0; i < *n; i++) {
	        if( !R_FINITE(NUMERIC_POINTER(sexp_fjac)[(*n)*j + i]) )
				error("Non-finite value(s) returned by jacobian (row=%d,col=%d)",i+1,j+1);
            rjac[(*ldr)*j + i] = NUMERIC_POINTER(sexp_fjac)[(*n)*j + i];
        }

    UNPROTECT(1);
}

SEXP nleqslv(SEXP xstart, SEXP fn, SEXP jac, SEXP rmethod, SEXP rglobal, SEXP rxscalm,
             SEXP control, SEXP rho)
{
    int     i, j;
    int     n, ldfjac, lr;
    int     info, nfev, njev;

    double  *x, *rwork, *rcdwrk, *qrwork;
	double  *xp, *fp, *gp, *scalex;
	int     *icdwrk, *outopt;
	const char *z;

    SEXP    eval_test;
    SEXP    sexp_x, sexp_diag, sexp_fvec, sexp_info, sexp_message, sexp_nfcnt, sexp_njcnt;
    SEXP    out, out_names;

    char    message[256];

	int     njcnt, nfcnt, termcd, lrwork, qrwsiz;
    int     maxit, jacflg, method, global, xscalm;
	double  xtol, ftol, btol, stepmx, dlt, sigma;

    PROTECT_INDEX ipx;

    if( activeflag )
        error("Recursive call of nleqslv not possible");
    ++activeflag;
        
    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));

    PROTECT(OS->x = duplicate(xstart));
    n = length(OS->x);

    switch (TYPEOF(OS->x)) {
    case REALSXP:
        break;
    case VECSXP:
        for (i = 0; i < n; i++)
            SET_VECTOR_ELT(OS->x, i, AS_NUMERIC(VECTOR_ELT(OS->x, i)));
        break;
    default:
        error("`x' that you provided is non-list and non-numeric!");
    }

    if (!isFunction(fn)) error("fn is not a function!");
    PROTECT(OS->fcall = lang2(fn, OS->x));

    if (!isEnvironment(rho)) error("rho is not an environment!");
    OS->env = rho;

    PROTECT(eval_test = eval(OS->fcall, OS->env));
    if (!IS_NUMERIC(eval_test))
        error("evaluation of fn function returns non-numeric vector!");
    i = length(eval_test);
    if( i != n )
        error("Length of fn result <> length of x!");

    for (i = 0; i < n; i++)
		if( !R_FINITE(NUMERIC_POINTER(eval_test)[i]) )
			error("evaluation of fn function has non-finite values\n   (starting at index=%d)",i+1);

    UNPROTECT(1);
    
    /*
     * query the optimal amount of memory Lapack needs
     * to execute blocked QR code
     * for largish n (500+) this speeds up significantly
     */
     
    F77_CALL(nwqmem)(&n, &qrwsiz);
    if( qrwsiz <= 0 )
        error("Error in querying amount of workspace for QR routines\n");
  
    qrwork  = real_vector(qrwsiz);
        
	lrwork  = 9*n+2*n*n;

    x       = real_vector(n);
    xp      = real_vector(n);
    fp      = real_vector(n);
    gp      = real_vector(n);
    scalex  = real_vector(n);
	rwork   = real_vector(lrwork);
	rcdwrk  = real_vector(3*n);
	icdwrk  = int_vector(n);
	outopt  = int_vector(2);

    xtol    = NUMERIC_VALUE(getListElement(control, "xtol"));
    ftol    = NUMERIC_VALUE(getListElement(control, "ftol"));
    btol    = NUMERIC_VALUE(getListElement(control, "btol"));
    sigma   = NUMERIC_VALUE(getListElement(control, "sigma"));
    stepmx  = NUMERIC_VALUE(getListElement(control, "stepmax"));
    dlt     = NUMERIC_VALUE(getListElement(control, "delta"));
    maxit   = INTEGER_VALUE(getListElement(control, "maxit"));

    outopt[0] = INTEGER_VALUE(getListElement(control, "trace"));
	outopt[1] = LOGICAL_VALUE(getListElement(control, "chkjac"));

	z = CHAR(STRING_ELT(rmethod, 0));
	if( strcmp(z,"Broyden") == 0 )
		method = 1;
	else
		method = 0;

	z = CHAR(STRING_ELT(rglobal, 0));
	if( strcmp(z,"qline") == 0 )
		global = 0;
	else if( strcmp(z,"gline") == 0 )
		global = 1;
	else if( strcmp(z,"dbldog") == 0 )
		global = 2;
	else if( strcmp(z,"pwldog") == 0 )
		global = 3;

	z = CHAR(STRING_ELT(rxscalm, 0));
	if( strcmp(z,"fixed") == 0 )
		xscalm = 0;
	else
		xscalm = 1;

    PROTECT_WITH_INDEX(sexp_diag = getListElement(control, "scalex"), &ipx);
    switch (TYPEOF(sexp_diag)) {
    case REALSXP:
        if (length(sexp_diag) == n) {
            REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
            for (i = 0; i < n; i++)
                scalex[i] = NUMERIC_POINTER(sexp_diag)[i];
        }
        else {
            REPROTECT(sexp_diag = NEW_NUMERIC(n), ipx);
        }
        break;
    case VECSXP:
        if (length(sexp_diag) == n) {
            REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
            for (i = 0; i < n; i++) {
                SET_VECTOR_ELT(sexp_diag, i, AS_NUMERIC(VECTOR_ELT(sexp_diag, i)));
                scalex[i] = NUMERIC_VALUE(VECTOR_ELT(sexp_diag, i));
            }
        }
        else {
            REPROTECT(sexp_diag = NEW_LIST(n), ipx);
            for (i = 0; i < n; i++)
                SET_VECTOR_ELT(sexp_diag, i, NEW_NUMERIC(1));
        }
        break;
    default:
        error("`scalex' that you provided is non-list and non-numeric!");
    }

    if (IS_NUMERIC(OS->x)) {
        for (i = 0; i < n; i++)
            x[i] = NUMERIC_POINTER(OS->x)[i];
    }
    else {
        for (i = 0; i < n; i++)
            x[i] = NUMERIC_VALUE(VECTOR_ELT(OS->x, i));
    }

/*========================================================================*/

    if( outopt[0] == 1)
        trace_header(method, global, xscalm, sigma, dlt, stepmx, ftol, xtol, btol);

    if (isNull(jac)) {
		jacflg = 0;
        F77_CALL(nwnleq)(x, &n, scalex, &maxit, &jacflg, &xtol, &ftol, &btol,
                         &method, &global, &xscalm, &stepmx, &dlt, &sigma,
                         rwork, &lrwork, rcdwrk, icdwrk, qrwork, &qrwsiz,
                         FCNJACDUM, &fcnval, outopt, xp, fp, gp, &njcnt, &nfcnt, &termcd);
    }
    else {
        if (!isFunction(jac))
            error("jac is not a function!");
        PROTECT(OS->jcall = lang2(jac, OS->x));
        PROTECT(eval_test = eval(OS->jcall, OS->env));
        if (!IS_NUMERIC(eval_test))
            error("evaluation of jac function returns non-numeric vector!");
        if (length(eval_test) != n*n)
            error("jac function must return numeric vector with length"
                  " == length(x) * length(x). Your function"
                  " returns one with length %d while %d expected.",
                  length(eval_test), n*n);
        UNPROTECT(1);

		jacflg = 1;
        F77_CALL(nwnleq)(x, &n, scalex, &maxit, &jacflg, &xtol, &ftol, &btol,
                         &method, &global, &xscalm, &stepmx, &dlt, &sigma,
                         rwork, &lrwork, rcdwrk, icdwrk, qrwork, &qrwsiz,
                         &fcnjac, &fcnval, outopt, xp, fp, gp, &njcnt, &nfcnt, &termcd);
        UNPROTECT(1);
    }

/*========================================================================*/

    fcn_message(message, termcd);

    PROTECT(sexp_x = NEW_NUMERIC(n));
    for (i = 0; i < n; i++)
        NUMERIC_POINTER(sexp_x)[i] = xp[i];

    PROTECT(sexp_fvec = NEW_NUMERIC(n));
    for (i = 0; i < n; i++)
        NUMERIC_POINTER(sexp_fvec)[i] = fp[i];

    PROTECT(sexp_info = NEW_INTEGER(1));
    INTEGER_POINTER(sexp_info)[0] = termcd;

    PROTECT(sexp_nfcnt = NEW_INTEGER(1));
    INTEGER_POINTER(sexp_nfcnt)[0] = nfcnt;

    PROTECT(sexp_njcnt = NEW_INTEGER(1));
    INTEGER_POINTER(sexp_njcnt)[0] = njcnt;

    PROTECT(sexp_message = NEW_STRING(1));
    SET_STRING_ELT(sexp_message, 0, mkChar(message));

    if (IS_NUMERIC(sexp_diag)) {
        for (i = 0; i < n; i++)
            NUMERIC_POINTER(sexp_diag)[i] = scalex[i];
    }
    else {
        for (i = 0; i < n; i++)
            NUMERIC_POINTER(VECTOR_ELT(sexp_diag, i))[0] = scalex[i];
    }

    PROTECT(out = NEW_LIST(7));
    SET_VECTOR_ELT(out, 0, sexp_x);
    SET_VECTOR_ELT(out, 1, sexp_fvec);
    SET_VECTOR_ELT(out, 2, sexp_info);
    SET_VECTOR_ELT(out, 3, sexp_message);
    SET_VECTOR_ELT(out, 4, sexp_diag);
    SET_VECTOR_ELT(out, 5, sexp_nfcnt);
    SET_VECTOR_ELT(out, 6, sexp_njcnt);

    PROTECT(out_names = NEW_STRING(7));
    SET_STRING_ELT(out_names, 0, mkChar("x"));
    SET_STRING_ELT(out_names, 1, mkChar("fvec"));
    SET_STRING_ELT(out_names, 2, mkChar("termcd"));
    SET_STRING_ELT(out_names, 3, mkChar("message"));
    SET_STRING_ELT(out_names, 4, mkChar("scalex"));
    SET_STRING_ELT(out_names, 5, mkChar("nfcnt"));
    SET_STRING_ELT(out_names, 6, mkChar("njcnt"));

    SET_NAMES(out, out_names);

    UNPROTECT(11);

    return out;
}

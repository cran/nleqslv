
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <R.h>

static	int jacsng = -1;
static  int jacupd = -1;
static  double jacond = 0.0;

/*
 * output for a single incorrect jacobian entry
 */
 
void F77_SUB(nwckot)(int *i, int *j, double *aij, double *wi)
{
    Rprintf("Chkjac  possible error in jacobian[%d,%d] = %20.13e\n"
            "                         Estimated[%d,%d] = %20.13e\n", *i, *j, *aij, *i, *j, *wi);
}

void F77_SUB(nwsnot)(int *jtype, int *ierr, double *rcond)
{
	/*
	 * save for later printing
	 */

	jacsng = *ierr;
	jacupd = *jtype;
	jacond = *rcond;
}

/*
 * output a compact description of the type of the jacobian used 
 *    Newton/Broyden followed by lowercase letter for ill-conditioned/singular
 *    estimated inverse condition number
 *
 * sprintf on Windows seems to use 3 digits for the exponent by default (msvcrt.dll?)
 * and doesn't obey the %7.1e format (uses 8.1)
 * that messes up a nice layout
 */

static void  nwrowhdr(int *iter)
{
	char jmethod;

    Rprintf( "  %4d ", *iter);
	if( jacupd < 0) { 
	    /* output padding */
        Rprintf("%11s","");
    }
    else {
    	jmethod = (jacupd == 0) ? 'N' : 'B';

        /* 
         * meaning jacsng
         *   0   jacobian is ok (not singular or ill-conditioned)
         *   1   jacobian is ill-conditioned
         *   2   jacobian is singular
         * 
         * Indicate this after output of <jmethod>
         */

    	if( jacsng == 0 )
    		Rprintf(" %c(%7.1e)", jmethod, jacond);
    	else if( jacsng == 1 )
    		Rprintf("%ci(%7.1e)", jmethod, jacond);
        else
    		Rprintf("%cs%9s", jmethod,"");

    	/*
    	 * avoid output of redundant information on next time called
    	 */

    	jacupd = -1;
	}
}

void F77_SUB(nwjerr)(int *iter)
{
    nwrowhdr(iter);
    Rprintf("\n");
}

static void enumout(double x)
{
    Rprintf(" %13.*e", fabs(x) >= 1E100? 5 : 6, x);
}

void F77_SUB(nwprot)(int *iter, int *lstep, double *oarg)
{
	double v;

	if( *lstep <= 0 ) {
		if( *lstep == -1)
			Rprintf("  %4s %11s %8s  %13s %13s\n",
						"Iter","Jac","Lambda","Fnorm","Largest |f|");

		Rprintf("  %4d%22s %13.6e %13.6e\n" , *iter, "", oarg[0],oarg[1]);
	}
	else {
		nwrowhdr(iter);
		v = *oarg;
		if( fabs(v) > 0.0001 )
			Rprintf( " %8.4f ",v);
		else
			Rprintf( " %8.1e ",v);

        enumout(oarg[1]);
        enumout(oarg[2]);
        Rprintf("\n");
	}
}

void F77_SUB(nwlsot)(int *iter, int *lstep, double *oarg)
{
	double v;

	if( *lstep <= 0 ) {
		if( *lstep == -1)
			Rprintf("  %4s %11s %8s  %13s %13s %13s\n",
						"Iter","Jac","Lambda","Ftarg","Fnorm","Largest |f|");

		Rprintf("  %4d%36s %13.6e %13.6e\n" , *iter, "", oarg[0],oarg[1]);
	}
	else {
		nwrowhdr(iter);
		v = *oarg;
		if( fabs(v) > 0.0001 )
			Rprintf( " %8.4f ",v);
		else
			Rprintf( " %8.1e ",v);

        enumout(oarg[1]);
        enumout(oarg[2]);
        enumout(oarg[3]);
        Rprintf("\n");
	}
}

void F77_SUB(nwdgot)(int *iter, int *lstep, double *oarg)
{
	char step;

	/*
 	 *  C gradient (cauchy) step
 	 *  N newton step
 	 *  P partial newton step
 	 *  W convex combination of P and C
     */

	if( *lstep <= 0 ) {
		if( *lstep == -1)
			Rprintf("  %4s %11s   %8s %8s %8s %8s %13s %13s\n",
 				       "Iter","Jac","Lambda", "Eta", "Dlt0", "Dltn", "Fnorm","Largest |f|");

            Rprintf("  %4d%50s" , *iter, "");
            enumout(oarg[0]);
            enumout(oarg[1]); 
            Rprintf("\n");
	}
	else {
		nwrowhdr(iter);
		step = "CWPN"[*lstep-1];
		Rprintf( " %c ", step);

		if( *lstep == 2 )
			Rprintf( "%8.4f", oarg[0]);
		else
			Rprintf( "%8s", "");

		Rprintf( " %8.4f %8.4f %8.4f",
                    oarg[3],oarg[1],oarg[2]);
        enumout(oarg[4]);
        enumout(oarg[5]);
        Rprintf("\n");
 	}
}

void F77_SUB(nwpwot)(int *iter, int *lstep, double *oarg)
{
	char step;

	/*
	 *  C gradient (cauchy) step
	 *  N newton step
	 *  W convex combination of P and C
	 */

	if( *lstep <= 0 ) {
		if( *lstep == -1)
			Rprintf("  %4s %11s   %8s %8s %8s %13s %13s\n",
						"Iter","Jac","Lambda", "Dlt0", "Dltn", "Fnorm","Largest |f|");

            Rprintf("  %4d%41s", *iter, "");
            enumout(oarg[0]);
            enumout(oarg[1]);
            Rprintf("\n");
	}
	else {
		nwrowhdr(iter);
		step = "CWN"[*lstep-1];
		Rprintf( " %c ", step);

		if( *lstep == 2 )
			Rprintf( "%8.4f",oarg[0]);
		else
			Rprintf( "%8s", "");

        Rprintf( " %8.4f %8.4f", oarg[1],oarg[2]);
        enumout(oarg[3]);
        enumout(oarg[4]);
        Rprintf("\n");
 	}
}

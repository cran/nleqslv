
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <R.h>

static	int jacsng = -1;
static  int jacupd = -1;
static  double jacond = 0.0;
static  double jacamu = 0.0;

void F77_SUB(nwstrot0)(const char *s, int *slen)
{   
    /*
     * called by Fortran subroutine to output <slen> chars from fortran char*(*)
     */
	int k, ns;

	ns = *slen;
	for(k=0; k < ns; k++)
		Rprintf("%c", s[k]);
	Rprintf("\n");
}

void F77_SUB(nwsnot)(int *jtype, int *ierr, double *rcond, double *mu)
{
	/*
	 * save for later printing
	 */

	jacsng = *ierr;
	jacupd = *jtype;
	jacond = *rcond;
	jacamu = *mu;
}

static char jcbuf[100];

static void  jackar()
{
	char jmethod;

	jcbuf[0] = 0;

	if( jacupd < 0) return;

	jmethod = (jacupd == 0) ? 'N' : 'B';

    /*
     * output inverse condition if jacobian is not singular
     * else output correction mu with an indication
     * of ill conditioning or singularity
     */

	if( jacsng == 0 )
		sprintf(jcbuf, " %c(%7.1e)", jmethod, jacond);
	else if( jacsng == 1 )
		sprintf(jcbuf, "%ci(%7.1e)", jmethod, jacamu);
    else
		sprintf(jcbuf, "%cs(%7.1e)", jmethod, jacamu);

	/*
	 * avoid output of redundant information on next time called
	 */

	jacupd = -1;

}

static void enumout(double x)
{
    Rprintf(" %13.*e", fabs(x) >= 1E100? 5 : 6, x);
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
		jackar();
		v = *oarg;
		Rprintf( "  %4d %11s", *iter, jcbuf);
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
			Rprintf("  %4s %11s   %8s %8s %8s %8s %8s %13s %13s\n",
 				       "Iter","Jac","Lambda","Gamma", "Eta", "Dlt0", "Dltn", "Fnorm","Largest |f|");

            Rprintf("  %4d%59s" , *iter, "");
            enumout(oarg[0]);
            enumout(oarg[1]); 
            Rprintf("\n");
	}
	else {
		jackar();
		step = "CNPW"[*lstep-1];
		Rprintf( "  %4d %11s %c ", *iter, jcbuf, step);

		if( *lstep == 4 )
			Rprintf( "%8.4f", oarg[0]);
		else
			Rprintf( "%8s", "");

		Rprintf( " %8.4f %8.4f %8.4f %8.4f",
                    oarg[3], oarg[4],oarg[1],oarg[2]);
        enumout(oarg[5]);
        enumout(oarg[6]);
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
		jackar();
		step = "CNW"[*lstep-1];
		Rprintf( "  %4d %11s %c ", *iter, jcbuf, step);

		if( *lstep == 3 )
			Rprintf( "%8.4f",oarg[0]);
		else
			Rprintf( "%8s", "");

        Rprintf( " %8.4f %8.4f", oarg[1],oarg[2]);
        enumout(oarg[3]);
        enumout(oarg[4]);
        Rprintf("\n");
 	}
}

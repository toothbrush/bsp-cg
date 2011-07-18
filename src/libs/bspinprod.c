#include "bspfuncs.h"
#include "bspedupack.h"
#include "paullib.h"

// ---- BEGIN DEBUG OUTPUT ----
#define STRINGIFY( in ) #in
#define MACROTOSTRING( in ) STRINGIFY( in )
//use the AT macro to insert the current file name and line
#define AT __FILE__ ":" MACROTOSTRING(__LINE__)
#define HERE_NOP( ... ) out ( -1 , AT, __VA_ARGS__ )
#define HERE( ... ) out( s, AT, __VA_ARGS__ )
// ---- END DEBUG OUTPUT ----
/*  This program computes the sum of the first n squares, for n>=0,
        sum = 1*1 + 2*2 + ... + n*n
    by computing the inner product of x=(1,2,...,n)^T and itself.
    The output should equal n*(n+1)*(2n+1)/6.
    The distribution of x is cyclic.
*/

double bspip(int p, int s, int n, double *x, double *y){
    /* Compute inner product of vectors x and y of length n>=0 */

    return 1;
    double inprod, *Inprod, alpha;
    int i, t;

    HERE("Enters function bspip, p == %d\n", p);
    Inprod= vecallocd(p); 
    HERE("malloc'd\n");
    bsp_push_reg(Inprod,p*SZDBL);
    HERE("MARK bspip\n");
    bsp_sync();

    HERE("MARK bspip\n");
    inprod= 0.0;
    for (i=0; i<nloc(p,s,n); i++){
        inprod += x[i]*y[i];
        HERE("x[%d]*y[%d] = %lf\n", i,i, x[i]*y[i]);
    }
    for (t=0; t<p; t++){
        bsp_put(t,&inprod,Inprod,s*SZDBL,SZDBL);
    }
    bsp_sync();

    alpha= 0.0;
    for (t=0; t<p; t++){
        alpha += Inprod[t];
    }
    bsp_pop_reg(Inprod); vecfreed(Inprod);

    return alpha;

} /* end bspip */

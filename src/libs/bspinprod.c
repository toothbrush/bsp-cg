#include "bspfuncs.h"
#include "bspedupack.h"

// This is my own version, bspip from BSPedupack
// cannot handle v1 and v2 having arbitrary distributions.

double bspip(int p,int s,int nv1, int nv2, double* v1, double *v2,
        int *procv2, int *indv2)
{
    /* Compute inner product of vectors x and y of length n>=0 */
    // here v1 is the local vector and v2 is the remote vector.

    double *v2_locals = vecallocd(nv1);

    int i;
    bsp_push_reg(v2, nv2*SZDBL);

    for(i=0; i<nv1; i++) {
        bsp_get(procv2[i], v2, indv2[i]*SZDBL, &v2_locals[i], SZDBL);

    }
    bsp_pop_reg(v2);

    double myip=0.0;

    for(i=0;i<nv1;i++) {
        myip += v1[i]*v2_locals[i];
    }

    double* Inprod = vecallocd(p);
    bsp_push_reg(Inprod, p*SZDBL);

    for(i=0;i<p;i++) {

        bsp_put(p, &myip, Inprod, s*SZDBL, SZDBL);

    }


    bsp_pop_reg(Inprod);

    double alpha = 0.0;
    for(i=0;i<p;i++)
        alpha += Inprod[i];

    free(Inprod);
    free(v2_locals);

    return alpha;

} /* end bspip */

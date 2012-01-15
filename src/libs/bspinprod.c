#include "bspfuncs.h"
#include "bspedupack.h"
#include "debug.h"

// This is my own version, bspip from BSPedupack
// cannot handle v1 and v2 having arbitrary distributions.

double bspip(int p,int s,int nv1, int nv2, double* v1, int*v1index,
                                                 double *v2,
                            int *procv2, int *indv2)
{
    /* Compute inner product of vectors x and y of length n>=0 */
    // here v1 is the local vector and v2 is the remote vector.

    double *v2_locals = vecallocd(nv1);

    int i;
    /*printf("registering address %p\n", v2);*/
    bsp_push_reg(v2, nv2*SZDBL);

    bsp_sync();
    for(i=0; i<nv1; i++) {
        //printf("bsp_get(procv2[i],\tv2,\tindv2[i]*SZDBL,\t&v2_locals[i],\tSZDBL);\n");
        //printf("bsp_get(%d,\t\tv2,\t%d*SZDBL,\t&v2_locals[%d],\tSZDBL);\n",procv2[i],indv2[i],i);
        bsp_get(procv2[v1index[i]], v2, indv2[v1index[i]]*SZDBL, &v2_locals[i], SZDBL);

        bsp_sync();
        // zero here is right if the call is from v.v... it's zero!
        //HERE("got %d = %lf\n", i, v2_locals[i]); // doesn't work if you don't sync first...
    }
    bsp_sync();
    bsp_pop_reg(v2);

    double myip=0.0;

    for(i=0;i<nv1;i++) {
        myip += v1[i]*v2_locals[i];
    }

    double* Inprod = vecallocd(p);
    bsp_push_reg(Inprod, p*SZDBL);
    bsp_sync();

    for(i=0;i<p;i++) {

        /*printf("bsp_put(i, &myip, Inprod, s*SZDBL, SZDBL);\n");*/
        /*printf("bsp_put(%d, &myip, Inprod, %d*SZDBL, SZDBL);\n",i,s);*/
        bsp_put(i, &myip, Inprod, s*SZDBL, SZDBL);

    }

    bsp_sync();

    bsp_pop_reg(Inprod);

    double alpha = 0.0;
    for(i=0;i<p;i++)
        alpha += Inprod[i];

    free(Inprod);
    free(v2_locals);

    return alpha;

} /* end bspip */

/*
 * copy distributed vec v into u
 */
void copyvec(int s,
        int nv, int nu,
        double* v, double* u,
        int* uindex,
        int* procv, int* indv)
{
    int i;

    bsp_push_reg(v, nv*SZDBL);
    bsp_sync();
    for(i=0;i<nu;i++) {

        HERE("bsp_get(procv[uindex[i]], v, indv[uindex[i]]*SZDBL, &u[i], SZDBL);\n");

        HERE("{\n\ti=%d\n\tuindex[i]=%d\n\tprocv=%d\n\tindv=%d\n}\n", i,uindex[i],procv[uindex[i]],indv[uindex[i]]);
        bsp_get(procv[uindex[i]], v, indv[uindex[i]]*SZDBL, &u[i], SZDBL);
        bsp_sync(); //useless
        HERE("we copied u[%d]=%lf\n",uindex[i],u[i]);

    }
    bsp_sync();
    bsp_pop_reg(v);

}

/*
 * add some other distributed vec to v (local)
 */

void addvec(int nv, double *v, int*vindex, int nr, double *remote,
        int *procr, int *indr) {

    double *tmp = vecallocd(nv);
    bsp_push_reg(remote,nr*SZDBL);
    bsp_sync();

    int i;
    for(i=0;i<nv;i++) {
        bsp_get(procr[vindex[i]], remote, indr[vindex[i]]*SZDBL, &tmp[i], SZDBL);
    }
    bsp_pop_reg(remote);
    bsp_sync();

    for(i=0;i<nv;i++) {
        v[i] += tmp[i];

    }

    free(tmp);
}

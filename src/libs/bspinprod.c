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

    double *v2_locals;
    v2_locals = vecallocd(nv1);

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
        HERE("got %d = %lf\n", i, v2_locals[i]); // doesn't work if you don't sync first...
    }
    bsp_sync();
    bsp_pop_reg(v2);

    double myip=0.0;

    for(i=0;i<nv1;i++) {
        myip += v1[i]*v2_locals[i];
    }

    double* Inprod;
    Inprod = vecallocd(p);
    HERE("Inprod = %p, p=%d\n", Inprod, p);
    bsp_push_reg(Inprod, p*SZDBL);
    bsp_sync();

    for(i=0;i<p;i++) {

        HERE("bsp_put(i, &myip, Inprod, s*SZDBL, SZDBL);\n");
        HERE("bsp_put(%d, &%f, Inprod, %d*SZDBL, SZDBL);\n",i,myip,s);

        bsp_put(i, &myip, Inprod, s*SZDBL, SZDBL);
        HERE("que?\n");

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
        int* vindex,
        int* procu, int* indu)
{
    int i;

    bsp_push_reg(u, nu*SZDBL);

    bsp_sync();


    for(i=0;i<nv;i++) {

        // put my v into u somewhere remote

        HERE("bsp_put(procu[vindex[i]], &v[i], u, indu[vindex[i]]*SZDBL, SZDBL);\n");
        HERE("{\n\ti=%d\n\tvindex[i]=%d\n\tprocu=%d\n\tindu=%d\n\tu=%p\n}\n",
                i,vindex[i],procu[vindex[i]],
                indu[vindex[i]],u);

        bsp_put(procu[vindex[i]], &v[i], u, indu[vindex[i]]*SZDBL, SZDBL);



    }


    bsp_sync();


    bsp_pop_reg(u);
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

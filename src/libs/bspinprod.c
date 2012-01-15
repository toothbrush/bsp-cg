#include "bspfuncs.h"
#include "bspedupack.h"
#include "debug.h"

// This is my own version; since bspip from BSPedupack
// cannot handle v1 and v2 having arbitrary distributions.

double bspip(int p,int s,
        int nv1, int nv2,
        double* v1, int*v1index,
        double *v2, int *procv2, int *indv2)
{
    /* Compute inner product of vectors x and y of length n>=0 */
    // here v1 is the local vector and v2 is the remote vector.

    double *v2_locals;
    v2_locals = vecallocd(nv1);

    int i;
    bsp_push_reg(v2, nv2*SZDBL);

    bsp_sync();
    for(i=0; i<nv1; i++) {
        // get all the vector components from v2, from where ever they're
        // stored. the array procv2 tells us this, but it is indexed by the
        // component's global index.
        bsp_get(procv2[v1index[i]], v2, indv2[v1index[i]]*SZDBL, &v2_locals[i], SZDBL);
    }

    double myip=0.0;

    bsp_sync();
    for(i=0;i<nv1;i++) {
        myip += v1[i]*v2_locals[i];
    }

#ifdef __GNUC__
    size_t tagsz;
#else
    int tagsz;
#endif
    tagsz = SZINT;
    bsp_set_tagsize(&tagsz);
    bsp_sync();

    for(i=0;i<p;i++) {

        if(s==i) // don't send myself messages.
            continue;

        bsp_send(i, &s, &myip, SZDBL);
    }

    bsp_sync();
    int nsums;
#ifdef __GNUC__
    size_t nbytes;
#else
    int nbytes;
#endif

    double alpha = myip;
    bsp_qsize(&nsums, &nbytes);
    int status, tag;
    bsp_get_tag(&status, &tag);

    for(i=0;i<nsums; i++) {
        bsp_move(&myip, SZDBL);

        alpha += myip;
        bsp_get_tag(&status, &tag);

    }

    bsp_pop_reg(v2);

    free(v2_locals);
    bsp_sync();

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
        // similar to how my bspip gets the proc that's relevant.
        bsp_get(procr[vindex[i]], remote, indr[vindex[i]]*SZDBL, &tmp[i], SZDBL);
    }
    bsp_pop_reg(remote);
    bsp_sync();

    for(i=0;i<nv;i++) {
        v[i] += tmp[i];

    }

    free(tmp);
}

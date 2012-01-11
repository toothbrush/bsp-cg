#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "libs/bspedupack.h"
#include "libs/bspfuncs.h"
#include "libs/vecio.h"
#include "libs/paullib.h"
#include "libs/debug.h"
// #include <Mondriaan.h>

#define EPS (10E-12)
#define KMAX (100)

/*
 * This program takes as input:
 *  - a matrix distributed over n processors
 *  - vector distributions u and v (which processor owns which value)
 *  - values for vector v (reals)
 *
 *  and outputs a vector u such that
 *  A . u == v
 *
 *  by using the conjugate gradient method, generalised to parallel
 *  as detailed in the documentation.
 */

int P;

char vfilename[STRLEN], ufilename[STRLEN], matrixfile[STRLEN];

void bspcg(){

    int s, p, n, nz, i, iglob, nrows, ncols, nv, nu,
        *ia, *ja, *rowindex, *colindex, *vindex, *uindex,
        *srcprocu, *srcindu, *destprocv, *destindv;
    double *a, *v, *u, *r, time0, time1, time2;

    bsp_begin(P);

    p= bsp_nprocs(); /* p=P */
    s= bsp_pid();

    // only proc 0 reads the files.
    if(s==0) {
        HERE("Start of BSP section.\n");
        char my_cwd[1024];
        getcwd(my_cwd, 1024);
        HERE("My working dir: PWD=%s\n", my_cwd);

        if(!file_exists(matrixfile)) {
            HERE("Matrix file doesn't exist. (%s)\n", matrixfile);
            bsp_abort(0);
        }
        if(!file_exists(vfilename)) {
            HERE("V-distrib file doesn't exist. (%s)\n", vfilename);
            bsp_abort(0);
        }
        if(!file_exists(ufilename)) {
            HERE("U-distrib file doesn't exist. (%s)\n", ufilename);
            bsp_abort(0);
        }
    }
    if (s==0){
        printf("CG solver\n");
        printf("   using %d processors\n",p);
    }

    /* Input of sparse matrix */
    bspinput2triple(matrixfile, p,s,&n,&nz,&ia,&ja,&a);
    HERE("Done reading matrix file.\n");

    /* Convert data structure to incremental compressed row storage */
    triple2icrs(n,nz,ia,ja,a,&nrows,&ncols,&rowindex,&colindex);
    HERE("Done converting to ICRS.\n");
    vecfreei(ja);

    /* Read vector distributions */
    bspinputvec(p,s,vfilename,&n,&nv,&vindex, &v);
    HERE("Loaded distribution vec v.\n");
    for(i=0; i<nv; i++){
        iglob= vindex[i];
    }

    bspinputvec(p,s,ufilename,&n,&nu,&uindex, &u);

    HERE("Loaded a %d*%d matrix, this proc has %d nz.\n", n,n,nz);
    if(s==0)
        printf("Loaded a %d*%d matrix, proc 0 has %d nz.\n", n,n,nz);

    assert(nu==nv); // we want the distribution to be equal

    if (s==0){
        HERE("Initialization for matrix-vector multiplications\n");
    }
    bsp_sync();
    time0= bsp_time();

    assert(ncols==nrows);
    // alloc metadata arrays
    srcprocu  = vecalloci(ncols);
    srcindu   = vecalloci(ncols);
    destprocv = vecalloci(nrows);
    destindv  = vecalloci(nrows);

    // do the heavy lifting.
    bsp_sync();
    time1= bsp_time();

    int k;

    k = 0; // iteration number
    r = vecallocd(nu);
    zero(nu,r);
    bspmv_init(p,s,n,nrows,ncols,nu,nv,rowindex,colindex,uindex,vindex,
               srcprocu,srcindu,destprocv,destindv);

    bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocu,srcindu,
            destprocv,destindv, nu, nv, u, r);
    negate(nv, r);
    axpy(nv, 1.0, v, r, r);

    long double rho = bspip(p,s,n,r,r);
    long double alpha,gamma,rho_old,beta;
    rho_old = 0; // just kills a warning.
    bsp_sync();

    double *pvec = vecallocd(nu);
    double *pold = vecallocd(nu);
    double *w    = vecallocd(nu);

    while ( k < KMAX &&
            rho > EPS * EPS * bspip(p,s,n,v,v)) {
        if ( k == 0 ) {
            copyvec(nv,r,pvec);
        } else {
            beta = rho/rho_old;
            axpy(nv,beta,pvec,r,     // beta*p + r
                              pvec); // into p
            if(s==0)
                printf("[Iteration %02d] rho  = %Le\n", k, rho);
        }
        bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocu,srcindu,
              destprocv,destindv,nu,nv,pvec,w);

        gamma = bspip(p,s,n,pvec,w);

        alpha = rho/gamma;

        axpy(nu,alpha,pvec,u,   // alpha*p + u
                           u);  // into u

        axpy(nv,-alpha,w,r,
                         r);

        rho_old = rho;
        rho     = bspip(p,s,n,r,r);

        k++;

    }

    // end heavy lifting.

    // postcondition:
    // u s.t. A.u = v

    bsp_sync();
    time2= bsp_time();

    if (s==0){
        HERE("End of matrix-vector multiplications.\n");
        HERE("Initialization took only %.6lf seconds.\n", time1-time0);
        printf("%d CG iterations took only %.6lf seconds.\n", k-1, (time2-time1));
        printf("The computed solution is:\n");
    }

    for(i=0; i<nu; i++){
        iglob=uindex[i];
        HERE("FINAL ANSWER *** proc=%d u[%d]=%lf \n",s,iglob,u[i]);
    }
    HERE("...which gives, filled in (should equal v):\n");
    bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocu,srcindu,
          destprocv,destindv,nu,nv,u,w);
    for(i=0; i<nu; i++){
        iglob=uindex[i];
        HERE("CHECKSUM     *** proc=%d A.u[%d]=%lf \n",s,iglob,w[i]);
    }

    double* answer = vecallocd(n);
    bsp_push_reg(answer,n*SZDBL);

    bsp_sync();

    for(i=0; i<nu; i++){
        iglob=uindex[i];
        bsp_put(0, &u[i], answer, iglob*SZDBL, SZDBL);
    }
    bsp_sync();

    if(s==0) {

        printf("========= Solution =========\n");
        printf("Final error = %Le\n\n", rho_old);

        for(i=0; i<n; i++) {
            printf("solution[%d] = %lf\n", i, answer[i]);
        }
    }

    bsp_pop_reg(answer);

    vecfreed(answer);
    vecfreed(w);        vecfreed(pvec);
    vecfreed(r);        vecfreed(pold);

    vecfreei(destindv); vecfreei(destprocv);
    vecfreei(srcindu);  vecfreei(srcprocu);
    vecfreed(u);        vecfreed(v);
    vecfreei(uindex);   vecfreei(vindex);
    vecfreei(rowindex); vecfreei(colindex);
    vecfreei(ia);       vecfreed(a);
    bsp_end();

} /* end bspcg */

int main(int argc, char **argv){

    bsp_init(bspcg, argc, argv);
    P = bsp_nprocs();

    if(argc != 3){
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t%s [mtx-dist] [v-vec]\n\n", argv[0]);
        exit(1);
    }

    strcpy(matrixfile, argv[1]);
    strcpy(vfilename, argv[2]);
    strcpy(ufilename, argv[2]);

    bspcg();
    exit(0);
}

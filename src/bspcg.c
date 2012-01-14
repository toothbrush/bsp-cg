#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "libs/bspedupack.h"
#include "libs/bspfuncs.h"
#include "libs/vecio.h"
#include "libs/paullib.h"
#include "libs/debug.h"

#define EPS (10E-12)
#define KMAX (1000)

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
        *ia, *ja, *rowindex, *colindex, *vindex, *uindex;
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
    HERE("Done converting to ICRS. nrows = %d, ncols = %d\n", nrows, ncols);
    assert(nrows == ncols); // since the matrix was square to start with.
    vecfreei(ja);

    /* Read vector distributions */
    bspinputvec(p,s,ufilename,&n,&nu,&uindex, &u);
    HERE("Loaded distribution vec u (nu=%d).\n",nu);
    for(i=0; i<nu; i++){
        iglob= uindex[i];
      //  HERE("original input vec %d = %lf\n", iglob, u[i]);
    }

    bspinputvec(p,s,vfilename,&n,&nv,&vindex, &v);
    HERE("Loaded distribution vec v (nv=%d). set to zero.\n",nv);
    zero(nv,v);

    HERE("Loaded a %d*%d matrix, this proc has %d nz.\n", n,n,nz);
    if(s==0)
        printf("Loaded a %d*%d matrix, proc 0 has %d nz.\n", n,n,nz);

    if (s==0){
        HERE("Initialization for matrix-vector multiplications\n");
    }
    bsp_sync();
    time0= bsp_time();

    // alloc metadata arrays
    int *srcprocv, *srcindv, *destprocu, *destindu;

    srcprocv  = vecalloci(ncols);
    srcindv   = vecalloci(ncols);
    destprocu = vecalloci(nrows);
    destindu  = vecalloci(nrows);

    // do the heavy lifting.
    bsp_sync();
    time1= bsp_time();

    int k;

    k = 0; // iteration number
    bspmv_init(p,s,n,nrows,ncols,nv,nu,rowindex,colindex,vindex,uindex,
               srcprocv,srcindv,destprocu,destindu);

    bsp_abort("normal...");

    r = vecallocd(nu);
    // corresponds to:
    // r := b - Ax,
    // but our guess for x = 0;
    for(i=0; i< nu; i++) {
        r[i] = u[i];
    }

    long double rho = bspip(p,s,nu,nu,r,r,destprocu,destindu);
    long double alpha,gamma,rho_old,beta;
    rho_old = 0; // just kills a warning.
    bsp_sync();

    double *pvec = vecallocd(nv);
    double *w    = vecallocd(nu);

    while ( k < KMAX &&
            rho > EPS * EPS * bspip(p,s,nv,nv,v,v,srcprocv,srcindv)) {
        if ( k == 0 ) {
            copyvec(nu, nv,r,pvec, srcprocv, srcindv);
        } else {
            beta = rho/rho_old;
            // TODO:
            axpy(nv,beta,pvec,r,     // beta*p + r
                              pvec); // into p
            if(s==0)
                printf("[Iteration %02d] rho  = %Le\n", k, rho);
        }
        bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
              destprocu,destindu,nv,nu,pvec,w);

        gamma = bspip(p,s,nv,nu,pvec,w,destprocu,destindu);

        alpha = rho/gamma;

        local_axpy(nu,alpha,pvec,u,   // alpha*p + u
                                 u);  // into u

        local_axpy(nu,-alpha,w,r,
                               r);

        rho_old = rho;
        rho = bspip(p,s,nu,nu,r,r,destprocu,destindu);

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
        printf("%d CG iterations took only %.6lf seconds (KMAX = %d).\n", k-1, (time2-time1), KMAX);
        printf("The computed solution is:\n");
    }

    for(i=0; i<nv; i++){
        iglob=vindex[i];
        HERE("FINAL ANSWER *** proc=%d v[%d]=%lf \n",s,iglob,v[i]);
    }
    HERE("...which gives, filled in (should equal v):\n");
    bspmv(p,s,n,nz,nrows,ncols,a,ia,destprocu,destindu,
          srcprocv,srcindv,nv,nu,v,w);
    for(i=0; i<nu; i++){
        iglob=uindex[i];
        HERE("CHECKSUM     *** proc=%d A.v[%d]=%lf \n",s,iglob,w[i]);
    }

    double* answer = vecallocd(n);
    bsp_push_reg(answer,n*SZDBL);
    int* nz_per_proc = vecalloci(P);
    bsp_push_reg(nz_per_proc,P*SZINT);

    bsp_sync();

    for(i=0; i<nv; i++){
        iglob=vindex[i];
        bsp_put(0, &v[i], answer, iglob*SZDBL, SZDBL);
    }
    bsp_put(0, &nz, nz_per_proc, s*SZINT, SZINT);
    bsp_sync();

    if(s==0) {

        int total_nz = 0;
        for(i=0; i<p;i++)
            total_nz += nz_per_proc[i];

        printf("========= Solution =========\n");
        printf("Final error = %Le\n\n", rho_old);
        printf("csv_answer:\tP\tN\tnz\ttime\titers\n");
        printf("csv_answer:\t%d\t%d\t%d\t%lf\t%d\n",P,n,total_nz,(time2-time1),k-1);

#ifdef DEBUG
        for(i=0; i<n; i++) {
            printf("solution[%d] = %lf\n", i, answer[i]);
        }
#endif
    }

    bsp_pop_reg(answer);
    bsp_pop_reg(nz_per_proc);

    vecfreed(answer);   vecfreei(nz_per_proc);
    vecfreed(w);        vecfreed(pvec);
    vecfreed(r);

    vecfreei(destindu); vecfreei(destprocu);
    vecfreei(srcindv);  vecfreei(srcprocv);
    vecfreed(u);        vecfreed(v);
    vecfreei(uindex);   vecfreei(vindex);
    vecfreei(rowindex); vecfreei(colindex);
    vecfreei(ia);       vecfreed(a);
    bsp_end();

} /* end bspcg */

int main(int argc, char **argv){

    bsp_init(bspcg, argc, argv);
    P = bsp_nprocs();

    if(argc != 4){
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t%s [mtx-dist] [u-dist] [v-dist]\n\n", argv[0]);
        exit(1);
    }

    strcpy(matrixfile, argv[1]);
    strcpy(ufilename, argv[2]);
    strcpy(vfilename, argv[3]);

    bspcg();
    exit(0);
}

#include <assert.h>
#include <string.h>
#include "libs/bspedupack.h"
#include "libs/bspfuncs.h"
#include "libs/vecio.h"
#include "libs/paullib.h"
#include "libs/debug.h"
// #include <Mondriaan.h>

#define EPS (10E-5)
#define KMAX (10)

#define DUMP( n, a ) for(counter0=0;counter0<n;counter0++) HERE("dump array[%d]=%lf\n",counter0, a[counter0])

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
        *srcprocv, *srcindv, *destprocu, *destindu;
    double *a, *v, *u, *r, time0, time1, time2;

    int counter0;

    bsp_begin(P);

    p= bsp_nprocs(); /* p=P */
    s= bsp_pid();

    HERE("Start of BSP section.\n");
    HERE("My working dir: PWD=%s\n", get_current_dir_name());

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
    /* Input of sparse matrix */
    bspinput2triple(matrixfile, p,s,&n,&nz,&ia,&ja,&a);
    HERE("Done reading matrix file.\n");
    // TODO: dump matrix A?
    HERE("First 5 values of a\n");
    DUMP(5,a);

    /* Convert data structure to incremental compressed row storage */
    triple2icrs(n,nz,ia,ja,a,&nrows,&ncols,&rowindex,&colindex);
    HERE("Done converting to ICRS.\n");
    vecfreei(ja);

    /* Read vector distributions */
    bspinputvec(p,s,vfilename,&n,&nv,&vindex, &v);
    HERE("Loaded distribution vec v. Dumping v =\n");
    DUMP(nv,v);
    bspinputvec(p,s,ufilename,&n,&nu,&uindex, &u);
    HERE("Loaded distribution vec u. Dumping u =\n");
    DUMP(nu,u);

    HERE("Some values: %d,%d,%d,%d,%d\n", p,s,n,nu,nv);

    assert(nu==nv); // right? we want the distribution to be equal?

    /*
    // this seems to work.
    // so is it true to say this program should work with all distributions?
    // i.e. r,p, etc should copy u/v's distr?
    for(i=0; i<nu; i++){
        iglob=uindex[i];
        HERE("proc=%d i=%d, u=%lf \n",s,iglob,u[i]);
    }
    for(i=0; i<nv; i++){
        iglob=vindex[i];
        HERE("proc=%d i=%d, v=%lf \n",s,iglob,v[i]);
    }
    bsp_abort(0);
    */

    if (s==0){
        HERE("CG solver\n");
        HERE(" using %d processors\n",p);
    }

    if (s==0){
        HERE("Initialization for matrix-vector multiplications\n");
    }
    bsp_sync();
    time0= bsp_time();

    // alloc metadata arrays
    srcprocv  = vecalloci(ncols);
    srcindv   = vecalloci(ncols);
    destprocu = vecalloci(nrows);
    destindu  = vecalloci(nrows);

    // do the heavy lifting.
    //EXAMPLE of Av: result goes into u.
    //bspmv_init(p,s,n,nrows,ncols,nv,nu,rowindex,colindex,vindex,uindex,
    //           srcprocv,srcindv,destprocu,destindu);
    //
    //bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
    //      destprocu,destindu,nv,nu,v,u);
    bsp_sync();
    time1= bsp_time();

    int k;

    k = 0; // iteration number
    r = vecallocd(nu);
    zero(nu,r);
    HERE("r. should be zero.\n");
    DUMP(nu,r);
    bspmv_init(p,s,n,nrows,ncols,nv,nu,rowindex,colindex,vindex,uindex,
               srcprocv,srcindv,destprocu,destindu);
    // must we not initialise?
    //
    // mv(A,u,r);
    // neg(r);
    // add(v,r,r);
    //
    // => r = v-A.u     ; the residue
    //

    // for a test: see if A.u gives something sensible.

    bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
            destprocu,destindu, nv, nu, u, r);
    HERE("Dumping A.u =\n");
    // TODO huh??? 0?
    DUMP(nu,r);
    negate(nv, r);
    axpy(nv, 1.0, v, r, r);

    HERE("Dumping v - A.u =\n");
    DUMP(nv,r);
    //for (i = 0 ; i < nv; i ++)
    //    HERE("uindex[%d]=%d\n", i, uindex[i]);

    printf("p %d, s %d, n %d, r %p\n", p,s,n,r);
    long double rho = bspip(p,s,n,r,r);
    long double alpha,gamma,rho_old,beta;
    rho_old = 0; // just kills a warning.
    bsp_sync();

    double *pvec = vecallocd(nu);
    double *pold = vecallocd(nu);
    double *w    = vecallocd(nu);

    while ( k < KMAX &&
            sqrt(rho) > EPS * sqrt(bspip(p,s,n,v,v))) {
        printf("sqrt(rho) = %lf\n", sqrt(rho));
        if ( k == 0 ) {
            copyvec(nv,r,pvec);
        } else {
            beta = rho/rho_old;
            axpy(nv,beta,pvec,r,     // beta*p + r
                              pvec); // into p
        }
        bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
              destprocu,destindu,nv,nu,pvec,w);

        gamma = bspip(p,s,n,pvec,w);

        alpha = rho/gamma;

        axpy(nu,alpha,pvec,u,   // alpha*p + u
                           u);  // into u

        axpy(nv,-alpha,w,r,
                         r);
        HERE("dumping R\n");
        DUMP(nv,r);

        rho_old = rho;
        rho     = bspip(p,s,n,r,r);

        k++;

        HERE("rho_old = %lf\n", rho_old);
        HERE("rho = %lf\n", rho);
        HERE("iteration %d\n", k);
    }

    // end heavy lifting.

    // postcondition:
    // u s.t. Au = v

    bsp_sync();
    time2= bsp_time();

    if (s==0){
        HERE("End of matrix-vector multiplications.\n");
        HERE("Initialization took only %.6lf seconds.\n",time1-time0);
        HERE("CG took only %.6lf seconds.\n",           (time2-time1));
        HERE("The computed solution is:\n");
    }

    for(i=0; i<nu; i++){
        iglob=uindex[i];
        HERE("FINAL ANSWER *** proc=%d u[%d]=%lf \n",s,iglob,u[i]);
    }
    HERE("...which gives, filled in (should equal v):\n");
    bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
          destprocu,destindu,nv,nu,u,w);
    for(i=0; i<nu; i++){
        iglob=uindex[i];
        HERE("FINAL ANSWER *** proc=%d A.u[%d]=%lf \n",s,iglob,w[i]);
    }


    vecfreed(w);        vecfreed(pvec);
    vecfreed(r);        vecfreed(pold);

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
        HERE_NOP("Usage:\n");
        HERE_NOP("\t%s [mtx-dist] [v-dist] [u-dist]\n\n", argv[0]);
        exit(1);
    }

    strcpy(matrixfile, argv[1]);
    strcpy(vfilename, argv[2]);
    strcpy(ufilename, argv[3]);

    bspcg();
    exit(0);
}

#include "libs/bspedupack.h"
#include "libs/bspfuncs.h"
#include <limits.h>
#include <stdlib.h>

#include "bspcg.h"
#include "paullib.h"

#define EPS (1.0E-12)
#define KMAX (100)

/* 
 * Author: Paul van der Walt
 * March 2011
 *
 * An implementation of parallel CG. 
 */

void bspInitCG(){
    
    double time0, time1;
    int p, s;
    
    bsp_begin(P);
    p= bsp_nprocs(); /* p = number of processors obtained */ 
    out("Now we have %d processors.\n", p);
    s= bsp_pid();    /* s = processor number */ 
    if (s==0){
        out("Hi.\n");
    }

    bsp_sync(); 
    time0=bsp_time();

    // begin work

    // input:  A: sparse n*n matrix
    //         b: dense vector, length n
    // output: x: dense vector, length n; A.x = b

    int n = 2 ; //TEMP TODO real $n$
    double *b;  // temp input vector

    // ICRS matrix A:
    int nz;     // number of nonzeroes
    int nrows;  // number of rows
    int ncols;  // number of columns
    double *a;  // the nonzero values
    int *inc;   // ??
    // end ICRS matrix A;
    // example data:
    nz = 4; nrows = 2; ncols = 2;
    a[0] = 4;    a[1] = 1;
    a[2] = 1;    a[3] = 3;
    inc[0] = 1;
    inc[0] = 1;
    inc[0] = 1 + ncols;
    inc[0] = 1;
    // end example data

    double *x = vecallocd(n);
    int i, k; 
    for(i = 0; i < n; i++)
    {
        x[i] = 0.0; // our best guess.
    }

    k = 0; // iteration number
    double* r = vecallocd(n);
    
    // temporary MV data structures
    int *srcprocv,*srcindv,*destprocu,*destindu;
    // end temporary MV ds

    bspmv_init(p,s,n, nrows, ncols, n,n, );
    r = bspmv(A,x);
    negate(r);
    add(r,r,b);
    double rho = bspip(p,s,n,r,r);
    double alpha,gamma,rho_old,beta;

    double *pvec = vecallocd(n);
    double *pold = vecallocd(n);
    double *w    = vecallocd(n);

    while(sqrt(rho) > EPS * sqrt(bspip(p,s,n,b,b)) && k < KMAX) {
        if(k == 0) {
            copyvec(pvec, r) ; // do p <- r;
        } else {
            beta = rho/rho_old;
            scalevec(beta,pold);
            add(pvec, r, pold); //TODO hmmmmmm p modified!
        }

        w = bspmv(A,pvec);
        gamma = bspip(p,s,n,pvec,w);
        alpha = rho / gamma;
        copyvec(pold, pvec);
        scalevector(pvec,alpha);
        scalevector(w,-alpha);
        add(x, x, pvec);
        add(r, r, w);
        rho_old = rho;
        rho = bspip(p,s,n,r,r);
        k++;
    }

    // here x should be correct

    // end work
    bsp_sync();  
    time1=bsp_time();

    fflush(stdout);
    fflush(stderr);
    if (s==0){
        out("This took only %.6lf seconds.\n", time1-time0);
        fflush(stdout);
        fflush(stderr);
    }

    bsp_end();

}

int main(int argc, char **argv){

    bsp_init(bspInitCG, argc, argv);

    /* sequential part */
    if (argc != 1)
    {
        out("Usage: %s\n", argv[0]);
        bsp_abort("This program expects an Extended Matrix Market format file on stdin. \n");
    }

    P = bsp_nprocs(); // maximum amount of procs

    out("Using %d processors. \n", P);

    /* SPMD part */
    bspInitCG();

    /* end, parallel part: sequential part */
    exit(0);

} /* end main */

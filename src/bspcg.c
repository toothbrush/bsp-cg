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

    double *x = vecallocd(n);
    int i, k; 
    for(i = 0; i < n; i++)
    {
        x[i] = 0.0; // our best guess.
    }

    k = 0; // iteration number
    double* r = vecallocd(n);
    r = axpy(A,x,minus(b));
    r = minus(r);
    double sqrtrho = norm(r); // the remaining residue
    double rho = sqrtrho*sqrtrho;

    while(sqrtrho > EPS * norm(b) && k < KMAX) {
        if(k == 0) {
            p = r;
        } else {
            beta = rho/rhoold;
            p = add(r,beta * p); //TODO axpy here.
        }

        w = axpy(A,p,0);
        gamma = vecmul(transpose(p),w);
        alpha = rho / gamma;
        x = add(x, alpha * p);
        r = add(r, -alpha * w);
        rhoold = rho;
        sqrtrho = norm(r);
        rho = sqrtrho * sqrtrho;
        k++;
    }

    // here x should be computed. 

    // end work
    bsp_sync();  
    time1=bsp_time();

    fflush(stdout);
    if (s==0){
        out("This took only %.6lf seconds.\n", time1-time0);
        fflush(stdout);
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

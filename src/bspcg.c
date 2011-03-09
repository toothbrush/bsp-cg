#include "bspedupack.h"
#include <limits.h>
#include <stdlib.h>

#include "bspcg.h"

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
    printf("Now we have %d processors.\n", p);
    s= bsp_pid();    /* s = processor number */ 
    if (s==0){
        printf("Hi.\n");
    }

    bsp_sync(); 
    time0=bsp_time();

    // begin work

    // end work
    bsp_sync();  
    time1=bsp_time();

    fflush(stdout);
    if (s==0){
        printf("This took only %.6lf seconds.\n", time1-time0);
        fflush(stdout);
    }

    bsp_end();

}

int main(int argc, char **argv){

    bsp_init(bspInitCG, argc, argv);

    /* sequential part */
    if (argc != 1)
    {
        printf("Usage: %s\n", argv[0]);
        bsp_abort("This program expects an Extended Matrix Market format file on stdin. \n");
    }

    P = bsp_nprocs(); // maximum amount of procs

    printf("Using %d processors. \n", P);

    /* SPMD part */
    bspInitCG();

    /* end, parallel part: sequential part */
    exit(0);

} /* end main */

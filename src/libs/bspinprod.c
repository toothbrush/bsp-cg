#include "bspfuncs.h"
#include "bspedupack.h"

/*  This program computes the sum of the first n squares, for n>=0,
        sum = 1*1 + 2*2 + ... + n*n
    by computing the inner product of x=(1,2,...,n)^T and itself.
    The output should equal n*(n+1)*(2n+1)/6.
    The distribution of x is cyclic.
*/

int P; /* number of processors requested */ 

double bspip(int p, int s, int n, double *x, double *y){
    /* Compute inner product of vectors x and y of length n>=0 */

    double inprod, *Inprod, alpha;
    int i, t;
  
    Inprod= vecallocd(p); bsp_push_reg(Inprod,p*SZDBL);
    bsp_sync();

    inprod= 0.0;
    for (i=0; i<nloc(p,s,n); i++){
        inprod += x[i]*y[i];
    }
    for (t=0; t<p; t++){
        bsp_put(t,&inprod,Inprod,s*SZDBL,SZDBL);
    }
    bsp_sync();

    alpha= 0.0;
    for (t=0; t<p; t++){
        alpha += Inprod[t];
    }
    bsp_pop_reg(Inprod); vecfreed(Inprod);

    return alpha;

} /* end bspip */

void bspinprod(){
    
    double *x, alpha, time0, time1;
    int p, s, n, nl, i, iglob;
    
    bsp_begin(P);
    p= bsp_nprocs(); /* p = number of processors obtained */ 
    s= bsp_pid();    /* s = processor number */ 
    if (s==0){
        out("Please enter n:\n"); fflush(stdout);
        scanf("%d",&n);
        if(n<0)
            bsp_abort("Error in input: n is negative");
    }
    bsp_push_reg(&n,SZINT);
    bsp_sync();

    bsp_get(0,&n,0,&n,SZINT);
    bsp_sync();
    bsp_pop_reg(&n);

    nl= nloc(p,s,n);
    x= vecallocd(nl);
    for (i=0; i<nl; i++){
        iglob= i*p+s;
        x[i]= iglob+1;
    }
    bsp_sync(); 
    time0=bsp_time();

    alpha= bspip(p,s,n,x,x);
    bsp_sync();  
    time1=bsp_time();

    out("Processor %d: sum of squares up to %d*%d is %.lf\n",
            s,n,n,alpha); fflush(stdout);
    if (s==0){
        out("This took only %.6lf seconds.\n", time1-time0);
        fflush(stdout);
    }

    vecfreed(x);
    bsp_end();

} /* end bspinprod */

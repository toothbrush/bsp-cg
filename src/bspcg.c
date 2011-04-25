#include <assert.h>
#include "libs/bspedupack.h"
#include "libs/bspfuncs.h"
#include "libs/vecio.h"
#include "libs/paullib.h"

#define EPS (1.0E-12)
#define KMAX (100)

/* This is a test program which uses bspmv to multiply a 
   sparse matrix A and a dense vector u to obtain a dense vector v.
   The sparse matrix and its distribution
   are read from an input file.
   The dense vector v is initialized by the test program.
   The distribution of v is read from an input file.
   The distribution of u is read from another input file.

   The output vector is defined by
       u[i]= (sum: 0<=j<n: a[i][j]*v[j]).
*/


int P;


void bspmv_test(){

    int s, p, n, nz, i, iglob, nrows, ncols, nv, nu, 
        *ia, *ja, *rowindex, *colindex, *vindex, *uindex,
        *srcprocv, *srcindv, *destprocu, *destindu;
    double *a, *v, *u, time0, time1, time2;
    char vfilename[STRLEN], ufilename[STRLEN], valuesfilename[STRLEN];

    bsp_begin(P);
    p= bsp_nprocs(); /* p=P */
    s= bsp_pid();


    
    /* Input of sparse matrix */
    bspinput2triple(p,s,&n,&nz,&ia,&ja,&a);

    /* Convert data structure to incremental compressed row storage */
    triple2icrs(n,nz,ia,ja,a,&nrows,&ncols,&rowindex,&colindex);
    vecfreei(ja);
    
    /* Read vector distributions */
    if (s==0){
        out("Please enter the filename of the v-vector distribution\n");
        scanf("%s",vfilename);
    }
    bspinputvec(p,s,vfilename,&n,&nv,&vindex);

    if (s==0){ 
        out("Please enter the filename of the u-vector distribution\n");
        scanf("%s",ufilename);
    }
    bspinputvec(p,s,ufilename,&n,&nu,&uindex);
    if (s==0){
        out("Sparse matrix-vector multiplication");
        out(" using %d processors\n",p);
    }

    /* Initialize input vector v */
    v= vecallocd(nv);

    /* Fill input vector with values */
    if (s==0){
        out("Finally enter the name of the file containing v's values\n");
        scanf("%s", valuesfilename);
    }
    readvalues(valuesfilename,nv,v);

    // postcondition: (not yet achieved)
    // u s.t. Au = v
    u= vecallocd(nu);
    
    if (s==0){
        out("Initialization for matrix-vector multiplications\n");
        fflush(stdout);
    }
    bsp_sync(); 
    time0= bsp_time();
    
    // alloc metadata arrays 
    srcprocv= vecalloci(ncols);
    srcindv= vecalloci(ncols);
    destprocu= vecalloci(nrows);
    destindu= vecalloci(nrows);

    assert(nu==nv);
    // do the heavy lifting.
    //EXAMPLE of Av: result goes into u.
    //bspmv_init(p,s,n,nrows,ncols,nv,nu,rowindex,colindex,vindex,uindex,
    //           srcprocv,srcindv,destprocu,destindu);

    //
    //EXAMPLE of Av: result goes into u. 
    //bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
    //      destprocu,destindu,nv,nu,v,u);
    bsp_sync(); 
    time1= bsp_time();

    int k;
    for(i = 0; i < nu; i++)
    {
        u[i] = 0.0; // our best guess.
    }

    k = 0; // iteration number
    double* r = vecallocd(nu);
    
    bspmv_init(p,s,n, nrows, ncols, nu,nu, rowindex,colindex,vindex,uindex, //input
           srcprocv, srcindv, destprocu, destindu ); // output
    bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,destprocu,destindu,nu,nu,u,r);
    negate(nu,r);
    addvec(nu,r,r,v);
    double rho = bspip(p,s,n,r,r);
    double alpha,gamma,rho_old,beta;

    double *pvec = vecallocd(n);
    double *pold = vecallocd(n);
    double *w = vecallocd(n);

    while(sqrt(rho) > EPS * sqrt(bspip(p,s,n,v,v)) && k < KMAX) {
        if(k == 0) {
            copyvec(nu, pvec, r) ; // do p <- r;
        } else {
            beta = rho/rho_old;
            scalevec(n,beta,pold);
            addvec(n,pvec, r, pold); //TODO hmmmmmm p modified!
        }
        out("Iteration %d.\n", k);

        bspmv_init(p,s,n, nrows, ncols, n,n, rowindex,colindex,vindex,uindex, //input
               srcprocv, srcindv, destprocu, destindu ); // output
        bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,destprocu,destindu,nu,nu,pvec,w);
        gamma = bspip(p,s,n,pvec,w);
        alpha = rho / gamma;
        copyvec(nu,pold, pvec);
        scalevec(n,alpha,pvec);
        scalevec(n,-alpha,w);
        addvec(n,u, u, pvec);
        addvec(n,r, r, w);
        rho_old = rho;
        rho = bspip(p,s,n,r,r);
        k++;
    }

    // here u should be correct
    
    // end heavy lifting.
    //
    bsp_sync();
    time2= bsp_time();
    
    if (s==0){
        out("End of matrix-vector multiplications.\n");
        out("Initialization took only %.6lf seconds.\n",time1-time0);
        out("CG took only %.6lf seconds.\n",           (time2-time1));
        out("The computed solution is:\n");
        fflush(stdout);
    }

    for(i=0; i<nu; i++){
        iglob=uindex[i];
        out("proc=%d i=%d, u=%lf \n",s,iglob,u[i]);
    }

    if (s==0) {
        out("Checksumming to follow.\n");
        out("A.u\t\tv (should be equal)\n");
    }

    bspmv_init(p,s,n, nrows, ncols, n,n, rowindex,colindex,vindex,uindex, //input
               srcprocv, srcindv, destprocu, destindu ); // output
    bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,destprocu,destindu,nu,nu,u,w);
    for(i=0; i<nu; i++) {
        out("%lf\t\t%lf\n",w[i],v[i]);
    }
    
    vecfreei(destindu); vecfreei(destprocu); 
    vecfreei(srcindv);  vecfreei(srcprocv); 
    vecfreed(u);        vecfreed(v);
    vecfreei(uindex);   vecfreei(vindex);
    vecfreei(rowindex); vecfreei(colindex); 
    vecfreei(ia);       vecfreed(a);
    bsp_end();
    
} /* end bspmv_test */

int main(int argc, char **argv){
 
    bsp_init(bspmv_test, argc, argv);
    P = bsp_nprocs();
    bspmv_test();
    exit(0);
}

#include "libs/bspedupack.h"
#include "libs/bspfuncs.h"
#include "libs/vecio.h"

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

#define NITERS 10

int P;


void bspmv_test(){

    int s, p, n, nz, i, iglob, nrows, ncols, nv, nu, iter,
        *ia, *ja, *rowindex, *colindex, *vindex, *uindex,
        *srcprocv, *srcindv, *destprocu, *destindu;
    double *a, *v, *u, time0, time1, time2;
    char vfilename[STRLEN], ufilename[STRLEN];

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
        printf("Please enter the filename of the v-vector distribution\n");
        scanf("%s",vfilename);
    }
    bspinputvec(p,s,vfilename,&n,&nv,&vindex);

    if (s==0){ 
        printf("Please enter the filename of the u-vector distribution\n");
        scanf("%s",ufilename);
    }
    bspinputvec(p,s,ufilename,&n,&nu,&uindex);
    if (s==0){
        printf("Sparse matrix-vector multiplication");
        printf(" using %d processors\n",p);
    }

    /* Initialize input vector v */
    v= vecallocd(nv);
    for(i=0; i<nv; i++){
        iglob= vindex[i];
        v[i]= iglob+1;
    }
    u= vecallocd(nu);
    
    if (s==0){
        printf("Initialization for matrix-vector multiplications\n");
        fflush(stdout);
    }
    bsp_sync(); 
    time0= bsp_time();
    
    srcprocv= vecalloci(ncols);
    srcindv= vecalloci(ncols);
    destprocu= vecalloci(nrows);
    destindu= vecalloci(nrows);
    bspmv_init(p,s,n,nrows,ncols,nv,nu,rowindex,colindex,vindex,uindex,
               srcprocv,srcindv,destprocu,destindu);

    if (s==0){
        printf("Start of %d matrix-vector multiplications.\n",(int)NITERS);
        fflush(stdout);
    }
    bsp_sync(); 
    time1= bsp_time();
    
    for(iter=0; iter<NITERS; iter++)
        bspmv(p,s,n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
              destprocu,destindu,nv,nu,v,u);
    bsp_sync();
    time2= bsp_time();
    
    if (s==0){
        printf("End of matrix-vector multiplications.\n");
        printf("Initialization took only %.6lf seconds.\n",time1-time0);
        printf("Each matvec took only %.6lf seconds.\n", 
                      (time2-time1)/(double)NITERS);
        printf("The computed solution is:\n");
        fflush(stdout);
    }

    for(i=0; i<nu; i++){
        iglob=uindex[i];
        printf("proc=%d i=%d, u=%lf \n",s,iglob,u[i]);
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
    printf("How many processors do you want to use?\n");
    scanf("%d",&P);
    if (P>bsp_nprocs()){
        printf("Not enough processors available:");
        printf(" %d wanted, %d available\n", P, bsp_nprocs());
        exit(1);
    }
    bspmv_test();
    exit(0);
}

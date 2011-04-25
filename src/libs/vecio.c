#include "vecio.h"
#include <stdio.h>
#include <limits.h>
#include "bspedupack.h"
#include "paullib.h"

// ---- BEGIN DEBUG OUTPUT ----
#define STRINGIFY( in ) #in
#define MACROTOSTRING( in ) STRINGIFY( in )
//use the AT macro to insert the current file name and line
#define AT __FILE__ ":" MACROTOSTRING(__LINE__)
#define HERE_NOP( ... ) out ( -1 , AT, __VA_ARGS__ )
#define HERE( ... ) out( s, AT, __VA_ARGS__ )
// ---- END DEBUG OUTPUT ----

void bspinput2triple(char*filename, int p, int s, int *pnA, int *pnz, 
                     int **pia, int **pja, double **pa){
  
    /* This function reads a sparse matrix in distributed
       Matrix Market format without the banner line
       from the input file and distributes
       matrix triples to the processors.
       The input consists of one line
           m n nz p  (number of rows, columns, nonzeros, processors)
       followed by p+1 lines with the starting numbers
       of the processor parts
           Pstart[0]
           Pstart[1]
           ...
           Pstart[p]
       which means that processor q will get all nonzeros
       numbered Pstart[q]..Pstart[q+1]-1.
       This is followed by nz lines in the format
           i j a     (row index, column index, numerical value).
       The input indices are assumed by Matrix Market to start
       counting at one, but they are converted to start from zero.
       The triples are stored into three arrays ia, ja, a,
       in arbitrary order.
       
       Input:
       p is the number of processors.
       s is the processor number, 0 <= s < p.

       Output:
       nA is the global matrix size.
       nz is the number of local nonzeros.
       a[k] is the numerical value of the k'th local nonzero,
            0 <= k < nz.
       ia[k] is the global row index of the  k'th local nonzero.
       ja[k] is the global column index.
    */

    int pA, mA, nA, nzA, nz, q, nzq, k, tagsz, status, *Pstart, *ia, *ja;
    double value, *a;
    indexpair t;
    FILE *fp;
    
    Pstart= vecalloci(p+1);
    bsp_push_reg(&nA,SZINT);
    bsp_push_reg(&nz,SZINT);
    tagsz= sizeof(indexpair);
    bsp_set_tagsize(&tagsz);
    bsp_sync();

    if (s==0){
        /* Open the matrix file and read the header */
        fp=fopen(filename,"r");

        /* A is an mA by nA matrix with nzA nonzeros
           distributed over pA processors. */
        fscanf(fp,"%d %d %d %d\n", &mA, &nA, &nzA, &pA);
        if(pA!=p)
            bsp_abort("Error: p not equal to p(A)\n"); 
        if(mA!=nA)
            bsp_abort("Error: matrix is not square");

        for (q=0; q<=p; q++)
            fscanf(fp,"%d\n", &Pstart[q]);
        for (q=0; q<p; q++){
            bsp_put(q,&nA,&nA,0,SZINT);
            nzq= Pstart[q+1]-Pstart[q];
            bsp_put(q,&nzq,&nz,0,SZINT);
        }
    }
    bsp_sync();


    /* Handle the processors one at a time.
       This saves buffer memory, at the expense of p-1 extra syncs.
       Buffer memory needed for communication is at most the maximum
       amount of memory a processor needs to store its vector components. */

    a= vecallocd(nz+1);
    ia= vecalloci(nz+1);  
    ja= vecalloci(nz+1);

    for (q=0; q<p; q++){      
        if (s==0){
            /* Read the nonzeros from the matrix file and
               send them to their destination */
            for (k=Pstart[q]; k<Pstart[q+1]; k++){
                fscanf(fp,"%d %d %lf\n", &t.i, &t.j, &value);
                /* Convert indices to range 0..n-1, 
                   assuming it was 1..n */
                t.i--;
                t.j--;
                /* Send a triple to P(q). Tag is a pair (i,j).
                   Payload is a numerical value */
                bsp_send(q,&t,&value,SZDBL);
            }
        }
        bsp_sync();
        
        if (s==q){
            /* Store the received nonzeros */
            for(k=0; k<nz; k++){
                bsp_get_tag(&status,&t);
                ia[k]= t.i;
                ja[k]= t.j;
                bsp_move(&a[k],SZDBL);
            }
        }
    }

    *pnA= nA;
    *pnz= nz;
    *pa= a;
    *pia= ia;
    *pja= ja;
    if (s==0)
        fclose(fp);
    bsp_pop_reg(&nz);
    bsp_pop_reg(&nA);
    bsp_sync();
    
} /* end bspinput2triple */
void triple2icrs(int n, int nz, int *ia,  int *ja, double *a,
                 int *pnrows, int *pncols,
                 int **prowindex, int **pcolindex){
    /* This function converts a sparse matrix A given in triple
       format with global indices into a sparse matrix in
       incremental compressed row storage (ICRS) format with 
       local indices.

       The conversion needs time and memory O(nz + sqrt(n))
       on each processor, which is O(nz(A)/p + n/p + p).
       
       Input:
       n is the global size of the matrix.
       nz is the local number of nonzeros.
       a[k] is the numerical value of the k'th nonzero
            of the sparse matrix A, 0 <= k <nz.
       ia[k] is the global row index of the k'th nonzero.
       ja[k] is the global column index of the k'th nonzero.
  
       Output:
       nrows is the number of local nonempty rows
       ncols is the number of local nonempty columns
       rowindex[i] is the global row index of the i'th
                   local row, 0 <= i < nrows.
       colindex[j] is the global column index of the j'th
                   local column, 0 <= j < ncols.
       a[k] is the numerical value of the k'th local nonzero of the
            sparse matrix A, 0 <= k < nz. The array is sorted by
            row index, ties being decided by column index.
       ia[k] = inc[k] is the increment in the local column index of the
              k'th local nonzero, compared to the column index of the
              (k-1)th nonzero, if this nonzero is in the same row;
              otherwise, ncols is added to the difference.
              By convention, the column index of the -1'th nonzero is 0.
   */
    
   int radix, i, iglob, iglob_last, j, jglob, jglob_last, k, inck,
       nrows, ncols, *rowindex, *colindex;
   
   /* radix is the smallest power of two >= sqrt(n)
      The div and mod operations are cheap for powers of two.
      A radix of about sqrt(n) minimizes memory and time. */

   for (radix=1; radix*radix<n; radix *= 2)
       ;
   
   /* Sort nonzeros by column index */
   sort(n,nz,ja,ia,a,radix,MOD);
   sort(n,nz,ja,ia,a,radix,DIV);
   
   /* Count the number of local columns */
   ncols= 0;
   jglob_last= -1;
   for(k=0; k<nz; k++){
       jglob= ja[k];
       if(jglob!=jglob_last)
           /* new column index */
           ncols++;
       jglob_last= jglob;
   }
   colindex= vecalloci(ncols);
   
   /* Convert global column indices to local ones.
      Initialize colindex */
   j= 0;
   jglob_last= -1;
   for(k=0; k<nz; k++){
       jglob= ja[k];
       if(jglob!=jglob_last){
           colindex[j]= jglob;
           j++;
       }
       ja[k]= j-1; /* local index of last registered column */
       jglob_last= jglob;
   }
   
   /* Sort nonzeros by row index using radix-sort */
   sort(n,nz,ia,ja,a,radix,MOD);
   sort(n,nz,ia,ja,a,radix,DIV);

   /* Count the number of local rows */
   nrows= 0;
   iglob_last= -1;
   for(k=0; k<nz; k++){
       iglob= ia[k];
       if(iglob!=iglob_last)
           /* new row index */
           nrows++;
       iglob_last= iglob;
   }
   rowindex= vecalloci(nrows);
                              
   /* Convert global row indices to local ones.
      Initialize rowindex and inc */
   i= 0;
   iglob_last= -1;
   for(k=0; k<nz; k++){
       if (k==0)
           inck= ja[k];
       else
           inck= ja[k] - ja[k-1];
       iglob= ia[k]; 
       if(iglob!=iglob_last){
           rowindex[i]= iglob;
           i++;
           if(k>0)
               inck += ncols;
       } 
       ia[k]= inck; /* ia is used to store inc */
       iglob_last= iglob;
   }
   if (nz==0)
       ia[nz]= 0;
   else 
       ia[nz]= ncols - ja[nz-1];
   ja[nz]= 0;                                                  
   a[nz]= 0.0;     
   
   *pncols= ncols;
   *pnrows= nrows;
   *prowindex= rowindex;
   *pcolindex= colindex;
   
} /* end triple2icrs */

void bspinputvec(int p, int s, const char *filename,
                 int *pn, int *pnv, int **pvindex,
                 double **pvalues){
  
    /* This function reads the distribution of a dense vector v
       from the input file and initializes the corresponding local
       index array.
       The input consists of one line
           n p    (number of components, processors)
       followed by n lines in the format
           i proc (index, processor number),
       where i=1,2,...,n.
       
       Input:
       p is the number of processors.
       s is the processor number, 0 <= s < p.
       filename is the name of the input file.

       Output:
       n is the global length of the vector.
       nv is the local length.
       vindex[i] is the global index corresponding to
                 the local index i, 0 <= i < nv.
    */

    int n, pv, q, np, b, i, k, globk, proc, ind, nv,
        *tmpproc, *tmpind, *Nv, *vindex;
    FILE *fp;

    double *allVals;

    bsp_push_reg(&n,SZINT);
    bsp_push_reg(&nv,SZINT);
    bsp_sync();

    if (s==0){
        /* Open the file and read the header */

        fp=fopen(filename,"r");
        fscanf(fp,"%d %d\n", &n, &pv);
        if(pv!=p)
            bsp_abort("Error: p not equal to p(vec)\n"); 
        for (q=0; q<p; q++)
            bsp_put(q,&n,&n,0,SZINT);
        allVals = vecallocd(n);  // place to temporarily store all vector reals. 
    }
    bsp_sync();

    /* The owner of the global index i and its local index
       are stored in temporary arrays which are distributed
       cyclically. */
    np= nloc(p,s,n);
    tmpproc= vecalloci(np);
    tmpind= vecalloci(np);
    bsp_push_reg(tmpproc,np*SZINT);
    bsp_push_reg(tmpind,np*SZINT);
    bsp_sync();

    if (s==0){
        /* Allocate component counters */
        Nv= vecalloci(p);
        for (q=0; q<p; q++)
            Nv[q]= 0;
    }

    /* block size for vector read */
    b= (n%p==0 ? n/p : n/p+1);
    for (q=0; q<p; q++){
        if(s==0){
            /* Read the vector components from file and
               put their owner and local index into their
               temporary location. This is done n/p components
               at a time to save memory  */
            for(k=q*b; k<(q+1)*b && k<n; k++){
                fscanf(fp,"%d %d %lf\n", &i, &proc, &allVals[k]); // also save value
                /* Convert index and processor number to ranges
                   0..n-1 and 0..p-1, assuming they were
                   1..n and 1..p */

                i--;  
                proc--;
                ind= Nv[proc];
                if(i!=k)
                    bsp_abort("Error: i not equal to index \n");
                bsp_put(i%p,&proc,tmpproc,(i/p)*SZINT,SZINT);
                bsp_put(i%p,&ind,tmpind,(i/p)*SZINT,SZINT);
                Nv[proc]++;
            }
        }
        bsp_sync();
    }

    if(s==0){
        for (q=0; q<p; q++)
            bsp_put(q,&Nv[q],&nv,0,SZINT);
        vecfreei(Nv);
    }
    bsp_sync();
    /* Store the components at their final destination */
    vindex= vecalloci(nv);  
    bsp_push_reg(vindex,nv*SZINT);
    bsp_sync();

    for(k=0; k<np; k++){
        globk= k*p+s;
        bsp_put(tmpproc[k],&globk,vindex,tmpind[k]*SZINT,SZINT);
    }
    bsp_sync();

    bsp_pop_reg(vindex);
    vecfreei(tmpind);
    vecfreei(tmpproc);

    /* grab reals from P0 */

    bsp_push_reg(allVals, n*SZDBL);
    bsp_sync();
    double *values = vecallocd(nv);

    for ( k = 0; k< nv; k++)
    {
        bsp_get(0, allVals, vindex[k]*SZDBL, &values[k], SZDBL);
    }
    bsp_pop_reg(allVals);
    bsp_sync();

    /* end value grabbing */

    *pn= n;
    *pnv= nv;
    *pvindex= vindex;
    *pvalues= values;

} /* end bspinputvec */

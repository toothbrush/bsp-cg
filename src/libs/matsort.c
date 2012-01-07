#include "matsort.h"
#include "vecio.h"
#include "bspedupack.h"

// This is from BSPedupack

int key(int i, int radix, int keytype){
   /* This function computes the key of an index i
      according to the keytype */
   
       if (keytype==DIV)
           return i/radix;
       else /* keytype=MOD */
           return i%radix;
           
} /* end key */



void sort(int n, int nz, int *ia, int *ja, double *a,
             int radix, int keytype){
   /* This function sorts the nonzero elements of an n by n sparse
      matrix A stored in triple format in arrays ia, ja, a.
      The sort is by counting. 
      If keytype=DIV, the triples are sorted by increasing value of
      ia[k] div radix.
      if keytype=MOD, the triples are sorted by increasing value of
      ia[k] mod radix.
      The sorting is stable: ties are decided so that the original
      precedences are maintained. For a complete sort by increasing
      index ia[k], this function should be called twice:
      first with keytype=MOD, then with keytype=DIV.
      
      Input: 
      n is the global size of the matrix.
      nz is the local number of nonzeros.
      a[k] is the numerical value of the k'th nonzero of the
           sparse matrix A, 0 <= k < nz.
      ia[k] is the global row index of the k'th nonzero.
      ja[k] is the global column index of the k'th nonzero.
      radix >= 1.
      
      Output: ia, ja, a in sorted order.
   */
   
   int *ia1, *ja1, nbins, *startbin, *lengthbin, r, k, newk;
   double *a1;
   
   ia1= vecalloci(nz);
   ja1= vecalloci(nz); 
   a1 = vecallocd(nz);
   
   /* Allocate bins */
   if (keytype==DIV)
       nbins= (n%radix==0 ? n/radix : n/radix+1);
   else if (keytype==MOD)
       nbins= radix;
   startbin= vecalloci(nbins);
   lengthbin= vecalloci(nbins);
       
   /* Count the elements in each bin */
   for (r=0; r<nbins; r++)
       lengthbin[r]= 0;
   for (k=0; k<nz; k++){
       r= key(ia[k],radix,keytype);
       lengthbin[r]++;
   }
    
   /* Compute the starting positions */
   startbin[0]= 0;
   for (r=1; r<nbins; r++)
       startbin[r]= startbin[r-1] + lengthbin[r-1];
       
   /* Enter the elements into the bins in temporary arrays (ia1,ja1,a1) */
   for (k=0; k<nz; k++){
       r= key(ia[k],radix,keytype);
       newk= startbin[r];
       ia1[newk]= ia[k];
       ja1[newk]= ja[k];
       a1[newk] = a[k];
       startbin[r]++;
   }
  
   /* Copy the elements back to the orginal arrays */
   for (k=0; k<nz; k++){
       ia[k]= ia1[k];
       ja[k]= ja1[k];
       a[k] = a1[k];
   }
   
   vecfreei(lengthbin);
   vecfreei(startbin);
   vecfreed(a1);
   vecfreei(ja1);
   vecfreei(ia1);
   
} /* end sort */

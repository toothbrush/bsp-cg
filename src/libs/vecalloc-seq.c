#include "vecalloc-seq.h"

/* These functions can be used to allocate and deallocate vectors and matrices.
   If not enough memory available, stop the program.

   This is a copy of the vecalloc?? functions from bspedupack.c,
   with the difference that bsp_abort isn't used here. Pity about the
   duplication, but oh well.
*/

double *vecallocd(int n){
    /* This function allocates a vector of doubles of length n */
    double *pd;

    if (n==0){
        pd= NULL;
    } else {
        pd= (double *)malloc(n*SZDBL);
        if (pd==NULL)
        {
            fprintf(stderr,"vecallocd: not enough memory\n");
            exit(20);
        }
    }
    return pd;

} /* end vecallocd */

ulong *vecalloculi(ulong n)
{
    /* This function allocates a vector of integers of length n */
    ulong *pi;

    if (n==0){
        pi= NULL; 
    } else { 
        pi= (ulong *)malloc(n*SZULL);
        if (pi==NULL)
        {
            fprintf(stderr,"vecalloculi: not enough memory\n");
            exit(20);
        }
    }
    return pi;

} /* end vecalloculi */
int *vecalloci(int n){
    /* This function allocates a vector of integers of length n */
    int *pi;

    if (n==0){
        pi= NULL; 
    } else { 
        pi= (int *)malloc(n*SZINT);
        if (pi==NULL)
        {
            fprintf(stderr,"vecalloci: not enough memory\n");
            exit(20);
        }
    }
    return pi;

} /* end vecalloci */

double **matallocd(int m, int n){
    /* This function allocates an m x n matrix of doubles */
    int i;
    double *pd, **ppd;

    if (m==0){
        ppd= NULL;  
    } else { 
        ppd= (double **)malloc(m*sizeof(double *));
        if (ppd==NULL)
        {
            fprintf(stderr,"matallocd: not enough memory\n");
            exit(20);
        }
        if (n==0){
            for (i=0; i<m; i++)
                ppd[i]= NULL;
        } else {  
            pd= (double *)malloc(m*n*SZDBL); 
            if (pd==NULL)
            {
                fprintf(stderr,"matallocd: not enough memory\n");
                exit(20);
            }
            ppd[0]=pd;
            for (i=1; i<m; i++)
                ppd[i]= ppd[i-1]+n;
        }
    }
    return ppd;

} /* end matallocd */

void vecfreed(double *pd){
    /* This function frees a vector of doubles */

    if (pd!=NULL)
        free(pd);

} /* end vecfreed */

void vecfreeuli(ulong *pi){
    /* This function frees a vector of integers */

    if (pi!=NULL)
        free(pi);

} /* end vecfreei */
void vecfreei(int *pi){
    /* This function frees a vector of integers */

    if (pi!=NULL)
        free(pi);

} /* end vecfreei */

void matfreed(double **ppd){
    /* This function frees a matrix of doubles */

    if (ppd!=NULL){
        if (ppd[0]!=NULL)
            free(ppd[0]);
        free(ppd);
    }

} /* end matfreed */

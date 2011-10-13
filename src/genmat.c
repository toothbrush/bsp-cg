#include "libs/vecalloc-seq.h"
#include <stdio.h>
#include <time.h>
#include "genmat.h"

#define true 1
#define false 0

typedef int bool;

int main (int argc, char** argv) {

    double sparsity;
    sparsity = 0.2; // nz = 20% of the size of the matrix

    int N;
    N = 10; // size of matrix.

    double mu;
    mu = 2.0; //factor for making matrix diagonal-dominant

    int nz = sparsity*N*N;
    int* xs;
    int* ys;
    double* vals;

    printf("Generating matrix. N=%d, density=%lf, target nz=%d\n",
            N, sparsity, nz);
    srand((unsigned)time(NULL));


    xs = vecalloci(nz);
    ys = vecalloci(nz);
    vals = vecallocd(nz);

    int v;


    for (v=0; v<nz; v++) {

        xs[v]= (double)N * ran();
        ys[v]= (double)N * ran();
        vals[v]= ran()*2-1; // [-1,1]

        printf("generated A[%d][%d]=%lf\n", xs[v],ys[v], vals[v]);

    }


    // oh my fucking god, look away, it's evil

    int uniques=0;
    bool found = false;

    int i;
    for(v=0;v<nz;v++) {
        // loop through the array, and for each element, check if the preceding 
        // portion of the array contains the same (i,j) pair. count first occurrences
        // of (i,j)'s.

        found=false;
        for(i=0;i<v;i++) {
            if(xs[v] == xs[i] &&
               ys[v] == ys[i]) {
                found = true;
                //break
            }
        }
        if(!found) {
            uniques ++;

        }

    }
    // now rebuild an array which actually contains those nonzeroes.

    int* fin_i;
    int* fin_j;
    double* fin_val;

    fin_i = vecalloci(uniques);
    fin_j = vecalloci(uniques);
    fin_val = vecallocd(uniques);
    uniques=0;
    for(v=0;v<nz;v++) {
        // loop through the array, and for each element, check if the preceding 
        // portion of the array contains the same (i,j) pair. count first occurrences
        // of (i,j)'s.

        found=false;
        for(i=0;i<v;i++) {
            if(xs[v] == xs[i] &&
               ys[v] == ys[i]) {
                found = true;
                //break
            }
        }
        if(!found) {
            // copy over.

            fin_i[uniques] = xs[v];
            fin_j[uniques] = ys[v];
            fin_val[uniques] = vals[v];
            printf("kept A[%d][%d]=%lf\n", fin_i[uniques],fin_j[uniques], fin_val[uniques]);

            uniques ++;

        }

    }






    printf("Generated %d nonzeroes\n", uniques);



    free(xs);
    free(ys);
    free(vals);
    free(fin_i);
    free(fin_j);
    free(fin_val);

    return 0;
}


double ran() {

    return ((double)rand()/(double)RAND_MAX);

}

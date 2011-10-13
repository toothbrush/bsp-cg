#include "libs/vecalloc-seq.h"
#include <stdio.h>


int main (int argc, char** argv) {

    double sparsity;
    sparsity = 0.2; // nz = 20% of the size of the matrix

    int N;
    N = 10; // size of matrix.

    double mu;
    mu = 2.0; //factor for making matrix diagonal-dominant

    int nz = sparsity*N;
    int* xs;
    int* ys;
    double* vals;

    printf("We'll now generate a matrix for you.\n");

    xs = vecalloci(nz);
    ys = vecalloci(nz);
    vals = vecallocd(nz);




    






    free(xs);
    free(ys);
    free(vals);

    return 0;
}

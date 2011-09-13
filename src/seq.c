#define N 2

#include "stdio.h"
#include "libs/bspedupack.h"

void cg_test();

int main(int argc, char** argv) {

    printf("hi\n");

    cg_test();

    return 0;

}





void cg_test() {

    int i;
    double **A;
    double *u;
    double *v;

    A = malloc(N * sizeof(double));

    u = malloc(N * sizeof(double));
    v = malloc(N * sizeof(double));

    for(i=0;i<N;i++)
    {
        A[i] = malloc(N * sizeof(double));
    }

    // fill A and v with values:

    A[0][0] = 4;
    A[0][1] = 1;
    A[1][0] = 1;
    A[1][1] = 3;




}


#define N 3
#define EPS (10E-8)
#define K_MAX (100)

#include "stdio.h"
#include "libs/bspedupack.h"
#include "libs/debug.h"

#define DUMP( n, a ) for(counter0=0;counter0<n;counter0++) printf("dump array[%d]=%lf\n",counter0, a[counter0])
void cg_test();

int main(int argc, char** argv) {

    printf("Starting CG, sequential.\n");

    cg_test();

    return 0;

}

void neg(double* x) {
    int i;
    for(i=0;i<N;i++)
        x[i] *= -1;
}

double ip(double* a, double* b) {

    int i;
    double tmp = 0;

    for(i=0;i<N;i++)
        tmp += a[i]*b[i];

    return tmp;
}

void scale(double factor, double* vec, double* res) {
    int i;

    for(i=0;i<N;i++)
        res[i] = factor * vec[i];

}

void add(double *a, double*b, double*res) {
    int i;

    for(i=0;i<N;i++)
        res[i] = a[i]+b[i];
}

void copy(double *in, double* out) {
    int i;

    for(i=0;i<N;i++)
        out[i] = in[i];
}

void mv(double** A, double *u, double*res) {

    int i,j;
    for(i=0;i<N;i++) {
        res[i] = 0;
        for(j=0;j<N;j++) {
            res[i] += A[i][j] * u[j];
            //printf("component A[%d][%d]*u[%d]=%lf\n", i,j,j,A[i][j] * u[j]);
        }
    }
}


void cg_test() {

    int i, counter0;
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

    A[0][0] = 1.0;
    A[0][1] = 2.0;
    A[0][2] = 3.0;
    A[1][0] = 4.0;
    A[1][1] = 5.0;
    A[1][2] = 6.0;
    A[2][0] = 7.0;
    A[2][1] = 8.0;
    A[2][2] = 9.0;

    v[0] = 11.0;
    v[1] = 12.0;
    v[2] = 13.0;

    //initialise guess of u:
    for(i=0;i<N;i++)
        u[i] = 0;

    // start 'heavy lifting'

    int k = 0; // iteration number

    double *r;
    r = malloc(N * sizeof(double));
    mv(A,u,r);
    neg(r);
    add(v,r,r);

    double rho, rho_old, beta, alpha, gamma;

    rho = ip(r,r);

    double *p; p = malloc(N*sizeof(double));
    double *w; w = malloc(N*sizeof(double));

    DUMP(N, r);
    while(sqrt(rho) > EPS * sqrt(ip(v,v)) &&
            k < K_MAX)
    {
        printf("sqrt(rho) = %lf\n", sqrt(rho));
        if(k==0) {
            copy(r,p);
        } else {
            beta = rho/rho_old;
            printf("is there a precision problem here? beta = rho/rho_old => %lf = %lf/%lf\n", beta, rho, rho_old);
            scale(beta, p, p);
            add(r,p,p);
        }
        mv(A,p,w);
        gamma = ip(p,w);
        alpha = rho/gamma;

        scale(alpha, p, p);
        add(u,p,u);

        neg(w);
        scale(alpha,w,w);
        add(r,w,r);
        rho_old = rho;
        rho = ip(r,r);
        k++;

        printf("iteration %d\n", k);

    }

    printf("my answer is:\n");
    for(i=0;i<N;i++)
        printf("u[%d] = %lf\n", i, u[i]);

    printf("Filling in gives:\n");
    mv(A,u,p);
    for(i=0;i<N;i++)
        printf("A.u[%d] = %lf\n", i, p[i]);

    free(p);
    free(w);
    free(r);
    free(u);
    free(v);
    for(i=0;i<N;i++)
        free(A[i]);
    free(A);


}


#define EPS (10E-8)
#define K_MAX (100)

#include "stdio.h"
#include "libs/bspedupack.h"
#include "libs/debug.h"

#define DUMP( n, a ) for(counter0=0;counter0<n;counter0++) printf("dump array[%d]=%lf\n",counter0, a[counter0])
void cg_test();
int N;

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

void find_N() {

    // temporary

    N = 10;
}

void read_A(double** A) {

    char* filename = "examples/10x10.mtx";
    int nz;

    FILE* fp = fopen(filename, "r");

    fscanf(fp, "%d %*d %d %*d", &N, &nz);

    printf("N = %d, nz = %d\n", N, nz);
    fscanf(fp, "%*d"); // ignore starting proc index
    fscanf(fp, "%*d"); // ignore ending   proc index

    int i, j; double val;
    int c;
    for(c=0;c<nz; c++) {
        fscanf(fp, "%d %d %lf", &i,&j,&val);
        A[i-1][j-1] = val;
    }

    fclose(fp);

}

void read_v(double* v) {

    char* filename = "examples/10x10.v";

    FILE* fp = fopen(filename, "r");

    fscanf(fp, "%*d %*d");

    int i; double val;
    int c;
    for(c=0;c<N; c++) {
        fscanf(fp, "%d %*d %lf", &i,&val);
        v[i-1] = val;
    }

    fclose(fp);


}

void cg_test() {

    int i,j, counter0=0;
    // remove warning:
    counter0=counter0;
    double **A;
    double *u;
    double *v;

    find_N();

    A = malloc(sizeof(double*)*N);
    for(i=0;i<N;i++)
        A[i] = malloc(sizeof(double)*N);

    // initialise A to all zeroes. (sparse)
    for(i = 0; i<N; i++)
        for(j = 0; j<N; j++)
            A[i][j] = 0.0;

    read_A(A);

    v = malloc(sizeof(double)*N);
    read_v(v);

    u = malloc(sizeof(double)*N);
    //initialise guess of u:
    for(i=0;i<N;i++)
        u[i] = 0.0;

    // start 'heavy lifting'

    int k = 0; // iteration number

    double *r;
    r = malloc(N * sizeof(double));
    mv(A,u,r);
    neg(r);
    add(v,r,r);

    double rho, rho_old, beta, alpha, gamma;

    rho_old = 0; //kill warning

    rho = ip(r,r);

    double *p; p = malloc(N*sizeof(double));
    double *w; w = malloc(N*sizeof(double));

    while(sqrt(rho) > EPS * sqrt(ip(v,v)) &&
            k < K_MAX)
    {
        if(k==0) {
            copy(r,p);
        } else {
            beta = rho/rho_old;
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

    }

    printf("after %d iterations, my answer is:\n", k);
    for(i=0;i<N;i++)
        printf("u[%d] = %lf\n", i, u[i]);

    printf("-> Filling in gives:\n");
    mv(A,u,p);
    for(i=0;i<N;i++)
        printf("A.u[%d] = %lf \t orig_v[%d]=%lf\n", i, p[i], i, v[i]);

    printf("\nFinal error=%lf\n", rho);
    free(p);
    free(w);
    free(r);
    free(u);
    free(v);
    for(i=0;i<N;i++)
        free(A[i]);
    free(A);

}

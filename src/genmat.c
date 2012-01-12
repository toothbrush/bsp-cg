#include <limits.h>
#include "libs/vecalloc-seq.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "genmat.h"
#include "libs/paulbool.h"

int N;

int main (int argc, char** argv) {

    double sparsity;
    // aim for a nonzero density given by sparsity:
    sparsity = 0.2; // nz = sparsity*100% of the size of the matrix

    /*
     * we say 'aim' here, since of course initially exactly
     *    nz = sparsity * N^2
     * nonzeroes will be generated at random spots, but because
     * the matrix must be symmetric and diagonally positive, the
     * actual number of nonzeroes will probably not be exactly
     * the projected number.
     */

    // read the desired size of the matrix from command line
    if (argc < 2) {
        printf("Usage: %s N [mu] [sparsity]\n", argv[0]);
        exit(-1);
    }

    if(sscanf(argv[1], "%d", &N) != 1) {
        printf("couldn't read command-line argument for N. must be an integer.\n");
        exit(-2);
    }
    double mu;
    mu = 2.5; //default scalar for making matrix diagonal-dominant

    // maybe the user supplied a different mu
    if(argc > 2 && sscanf(argv[2], "%lf", &mu) != 1) {
        exit(-2);
    }
    // maybe the user supplied a different sparsity
    if(argc > 3 && sscanf(argv[3], "%lf", &sparsity) != 1) {
        exit(-2);
    }

    int nz = sparsity*N*N;
    int* xs;
    int* ys;
    double* vals;

    fprintf(stderr,"Generating matrix. N=%d, density=%lf, target nz=%d, ",
            N, sparsity, nz);
    fprintf(stderr, "mu = %lf\n", mu);

    // seed the random generator.
    srandom((unsigned)time(NULL));

    fprintf(stderr,"INT_MAX = %d\n", INT_MAX);

    xs = vecalloci(nz);
    ys = vecalloci(nz);
    vals = vecallocd(nz);

    bool* diag_done;
    diag_done = malloc(N*sizeof(bool));

    int i;
    for(i = 0; i<N; i++) {
        diag_done[i] = false;
    }

    int nz_generated;

    int x,y;
    nz_generated = 0;
    double fake_transpose;
    x=0;y=0;
    while(x<N && y<N) { //don't escape matrix bounds.

        fprintf(stderr,"progress: %lf%%\n", 100*nz_generated/(double)(nz/2.0));
        if(x==y) {
            //diagonal, so always generate.

            xs[nz_generated]=x;
            ys[nz_generated]=y;
            vals[nz_generated]=ran()*2.0-1.0;
            diag_done[x] = true;

#ifdef DEBUG
            fprintf(stderr,"generated A[//][%d]=%lf\n"
                                                      , ys[nz_generated]
                                                      , vals[nz_generated]);
#endif

            nz_generated++;

        } else {
            // not a diagonal. only add if in
            // lower triangular half
            if(x<y) {

                if(nz_generated > nz) {
                    // this should NEVER happen, although all that's
                    // stopping it from happening is our ran() being
                    // well-behaved...........
                    printf("EEK! something went wrong!!\n");
                    exit(666);
                }

                xs[nz_generated]= x;
                ys[nz_generated]= y;
                // simulate the distribution of values which
                // would occur if we do A+A^T afterwards.
                if (ran() < sparsity) {
                    fake_transpose    = ran()*2.0-1.0;
                    vals[nz_generated]= ran()*2.0-1.0 + fake_transpose;
                } else {
                    vals[nz_generated]= ran()*2.0-1.0;
                }
#ifdef DEBUG
            fprintf(stderr,"generated A[%d][%d]=%lf\n", xs[nz_generated]
                                                      , ys[nz_generated]
                                                      , vals[nz_generated]);
#endif
                nz_generated++;


            }
        }
        x += 1/sparsity * (ran() + 0.5);
        if( x >= N ) {
            y += x/N;
            x  = x%N;
        }

    }

    fprintf(stderr, "generated initial randoms\n");

    int diagonals_present = 0;
    for(i=0; i<nz_generated; i++) {
        if(xs[i]==ys[i])
            diagonals_present++;
    }

    fprintf(stderr,"generated %d nzeros, array was %d big.\n", nz_generated, nz);


#ifdef DEBUG
    fprintf(stderr, "found %d diagonal(s), still need %d more.\n", diagonals_present, (N-diagonals_present));
#endif

    // add the missing diagonals, and add mu to each diagonal.
    int newsize = nz_generated + (N - diagonals_present);
    fprintf(stderr,"reallocating arrays to %d long\n", newsize);
    int* diag_i;
    int* diag_j;
    double* diag_val;

    diag_i   = realloc(xs  ,SZINT*newsize);
    diag_j   = realloc(ys  ,SZINT*newsize);
    diag_val = realloc(vals,SZDBL*newsize);

    if(diag_i == NULL ||
            diag_j == NULL ||
            diag_val == NULL)
    {
        printf("out of memory!");
        exit(44);
    }

    fprintf(stderr, "newsize = %d\n", newsize);
    addDiagonal(mu, diag_i, diag_j, diag_val, nz_generated, newsize, diag_done);
    nz_generated=newsize;
#ifdef DEBUG
    for(i=0;i<newsize;i++) {
        fprintf(stderr,"after addDiagonal A[%d][%d]=%lf\n", diag_i[i],diag_j[i], diag_val[i]);
    }

    fprintf(stderr, "Going to make symmetric now... (nz_generated = %d)\n", nz_generated);
#endif

    // now we explicitly fill the array with the
    // upper triangle values

    // things must be symmetric, but they aren't, yet
    // ... here's a good place to do the transposing thing.

    newsize = nz_generated * 2 - N; //number of real nonzeros, don't
                                    // count diagonals twice.
    int *new_i;
    int *new_j;
    double *new_v;

    fprintf(stderr, "newsize = %d\n", newsize);
    fprintf(stderr, "%p\n", diag_j);
    new_i = realloc(diag_i  ,SZINT*newsize);
    fprintf(stderr, "new_j = realloc(diag_j  ,SZINT*newsize);\n");
    new_j = realloc(diag_j  ,SZINT*newsize);
    fprintf(stderr, "%p = realloc(%p  ,%lu*%d);\n", new_j, diag_j, SZINT,newsize);
    new_v = realloc(diag_val,SZDBL*newsize);
    fprintf(stderr, "new i and j = %p, %p, %p\n", new_i, new_j, new_v);

    if(new_i == NULL ||
            new_i == NULL ||
            new_v == NULL)
    {
        printf("out of memory (2)!");
        exit(44);
    }
    diag_i = new_i;
    diag_j = new_j;
    diag_val = new_v;
    addTranspose(newsize,diag_i,diag_j,diag_val,
                              nz_generated);

#ifdef DEBUG
    for(i=0;i<newsize;i++)
        // to make diags stand out.
        if(diag_i[i]==diag_j[i])
            fprintf(stderr,"after transpose A[%d][%d]=%lf \\\\\n", diag_i[i],diag_j[i], diag_val[i]);
        else
            fprintf(stderr,"after transpose A[%d][%d]=%lf\n", diag_i[i],diag_j[i], diag_val[i]);
#endif

    checkStrictDiagonallyDominant(diag_i,diag_j,diag_val, newsize);

    // now quickly generate a test-vector to solve against:

    double *vec = vecallocd(N);
    for(i=0;i<N;i++)
        vec[i]=ran();

    fprintf(stderr,"Left with %d nonzeroes; nonzero density = %lf\n", newsize, newsize/((double)N*N));
    fprintf(stderr,"========== OUTPUTTING ... ==========\n");

    outputSimpleMatrix(newsize, diag_i, diag_j, diag_val, vec);
    //outputMatrix(newsize, diag_i, diag_j, diag_val, vec);
    //outputMathematicaMatrix(newsize, diag_i, diag_j, diag_val, vec);

    free(diag_done);
    free(vec);
    free(diag_i);
    free(diag_j);
    free(diag_val);

    return 0;
}

void addDiagonal(double mu, int* i, int* j, double* v, int nz_generated, int newsize, bool* diags_done) {

    int c;
    int pos = nz_generated;
    for(c=0; c<N; c++) {
        if(!diags_done[c]){

            i[pos] = c;
            j[pos] = c;
            v[pos] = ran()*2.0-1.0;

            fprintf(stderr, "adding diag %d to array pos %d\n", c, pos);

            pos++;


        }

    }
    // now add mu to each diagonal value.

    fprintf(stderr, "doing for c=0;c<%d;c++\n", newsize);

    for(c=0;c<newsize; c++) {

        if(i[c]==j[c]) {
            v[c] += mu;

        }

    }


}

/**
 * This function takes some matrix A and produces (almost) A' = A + A^T, which is
 * a symmetric matrix. Note that the diagonals are not added.
 *
 */
void addTranspose(int final_size, int* i, int* j, double* v,
                 int nz_generated) {

    int pos = nz_generated;

    int c;
    for(c=0;c<nz_generated;c++) {

        if(i[c] != j[c]) {
            // not a diagonal, so copy over.

            i[pos] = j[c];
            j[pos] = i[c];
            v[pos] = v[c];

            pos++;

        }
    }

}


void checkStrictDiagonallyDominant(int* i, int* j, double* v, int nz)
{

    // steps:
    // first sum all rows
    // then find diagonals
    // check each diagonal against the summed rows.

    int c;
    double * rowtotal;
    rowtotal = vecallocd(N);
    double * diagonals;
    diagonals = vecallocd(N);

    for(c = 0; c< N; c++)
    {
        rowtotal[c] = 0;
        diagonals[c] = 0;
    }

    for(c = 0; c< nz; c++)
    {
        // find diagonals:
        if(i[c] == j[c]){
            diagonals[i[c]] = v[c];
        } else {
            rowtotal[i[c]] += fabs(v[c]);
        }

    }

    // foreach diag, check.
    for(c=0; c<N; c++) {

        if ( !(fabs(diagonals[c]) > rowtotal[c]) ) {
            fprintf(stderr, "PROBLEM: diagonal > rowtotal doesn't hold: \n"
                            "    diagonals[%d] = %lf\n"
                            "    rowtotal[%d]  = %lf\n",
                            c, fabs(diagonals[c]),
                            c, rowtotal[c]
                   );
            fprintf(stderr, "increase mu? sometimes just running again is enough.\n");
            exit(5);
        }

    }

    free(rowtotal);
    free(diagonals);

}

/*
 * return a random double in the interval [0,1]
 */
double ran() {

    return ((double)random()/(double)RAND_MAX);

}

void outputMathematicaMatrix(int nz, int*i, int*j, double*v, double*vec) {

    // here we'll print the matrix in Mathematica format

    // Mathematica "header"

    printf("somemat = SparseArray[ { \n");
    // the value lines: i j value:
    int c;
    for(c=0;c<nz-1;c++) {

        printf("{%d, %d} -> %lf,\n", i[c]+1, j[c]+1, v[c]); // Mathematica expects 1-based coordinates.

    }
    printf("{%d, %d} -> %lf\n", i[nz-1]+1, j[nz-1]+1, v[nz-1]); // Mathematica expects 1-based coordinates.

    printf("} ] ;\n");
    printf("(* somemat // MatrixForm *)\n");

    printf("\n\n\n");

    printf("(* ======= vector v follows ====== *)\n");


    printf("vec = {\n");
    // N vector entries, one proc.
    for(c=0;c<N-1;c++) {
        // each line is a value, in order of the vector indices.
        printf("%lf,\n", vec[c]);
    }
    // last line without comma.
    printf("%lf\n};\n(* vec // MatrixForm *)\n", vec[N-1]);

    // and finally, for the paranoid:
    printf("\n\nPositiveDefiniteMatrixQ[somemat]\n\nSymmetricMatrixQ[somemat]\n");
    printf("correctAnswer = LinearSolve[somemat,vec];\n");
    printf("(* correctAnswer // MatrixForm *)\n");

}
void outputSimpleMatrix(int nz, int*i, int*j, double*v, double*vec) {

    // here we'll print the matrix in simple format, for one proc.
    // this is nice for testing CG without going through Mondriaan first.

    // no header.
    //size line: m n nz
    printf("%d %d %d %d\n", N, N, nz, 1); // one processor, i.e. not distributed
    // begin and end bounds of proc 1 are 0 and nz
    printf("0\n%d\n", nz);

    // next the value lines: i j value:
    int c;
    for(c=0;c<nz;c++) {

        printf("%d %d %lf\n", i[c]+1, j[c]+1, v[c]); // bsp-cg expects 1-based coordinates.

    }

    fprintf(stderr, "======= vector v follows ======\n");


    // N vector entries, one proc.
    printf("%d %d\n", N, 1);
    for(c=0;c<N;c++) {
        // each line is
        //   coordinate,processor,value
        printf("%d %d %lf\n", c+1, 1, vec[c]);
    }


}
void outputMatrix(int nz, int*i, int*j, double*v, double*vec) {

    // here we'll print the matrix in EMM format.

    //header:
    printf("%%%%Extended-MatrixMarket matrix coordinate double general original\n");
    //size line: m n nz
    printf("%d %d %d\n", N, N, nz);

    // next the value lines: i j v,
    int c;
    for(c=0;c<nz;c++) {

        printf("%d %d %lf\n", i[c]+1, j[c]+1, v[c]); // Mondriaan expects 1-based coordinates.

    }

    // ...and now for the vector-to-be-solved-for:
    printf("%%%%b vector double general array original\n");
    //size line:
    printf("%d\n", N);
    for(c=0;c<N;c++) {
        printf("%lf\n", vec[c]);
    }


}


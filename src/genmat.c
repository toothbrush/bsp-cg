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


    xs = vecalloci(nz);
    ys = vecalloci(nz);
    vals = vecallocd(nz);

    int v,i,nz_generated;

    int x,y;
    nz_generated = 0;
    double fake_transpose;
    // try each matrix location and generate a value with some
    // probability.
    for(i=0; i< N*N; i++) {

        x=i/N;
        y=i%N;

        if(x==y) {
            //diagonal, so always generate.
            //nice, guarantee that all diagonals are set.

            xs[nz_generated]=x;
            ys[nz_generated]=y;
            vals[nz_generated]=ran()*2.0-1.0;

            // and it's a diagonal, so add MU
            vals[nz_generated] += mu;
#ifdef DEBUG
            fprintf(stderr,"generated A[//][%d]=%lf\n"
                                                      , ys[nz_generated]
                                                      , vals[nz_generated]);
#endif

            nz_generated++;

        } else {
            // not a diagonal. only add if in
            // lower triangular half
            if(x<y && sparsity > ran()) {

                xs[nz_generated]= i/N;
                ys[nz_generated]= i%N;
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


    }

    fprintf(stderr, "generated initial randoms\n");

    // oh my fucking god, look away, it's evil

    int diagonals_present = 0;
    for(i=0; i<nz_generated; i++) {
        if(xs[i]==ys[i])
            diagonals_present++;
    }

    printf("generated %d nzeros, array was %d big.\n", nz_generated, nz);


#ifdef DEBUG
    fprintf(stderr, "found %d diagonal(s), still need %d more.\n", diagonals_present, (N-diagonals_present));
#endif

    // now we explicitly fill the array with the
    // upper triangle values

    int newsize = nz_generated * 2 - N; //number of real nonzeros, don't
                                        // count diagonals twice.
    int* diag_i;
    int* diag_j;
    double* diag_val;

    diag_i = realloc(xs,SZINT*newsize);
    diag_j = realloc(ys,SZINT*newsize);
    diag_val = realloc(vals,SZDBL*newsize);

    if(diag_i == NULL ||
            diag_j == NULL ||
            diag_val == NULL)
    {
        printf("out of memory!");
        exit(44);
    }

#ifdef DEBUG
    fprintf(stderr, "Going to make symmetric now...\n");
#endif

    // things must be symmetric, but they aren't, yet
    // ... here's a good place to do the transposing thing.


    addTranspose(newsize,diag_i,diag_j,diag_val,
                              nz_generated);

#ifdef DEBUG
    for(v=0;v<newsize;v++)
        // to make diags stand out.
        if(diag_i[v]==diag_j[v])
            fprintf(stderr,"after transpose A[%d][%d]=%lf \\\\\n", diag_i[v],diag_j[v], diag_val[v]);
        else
            fprintf(stderr,"after transpose A[%d][%d]=%lf\n", diag_i[v],diag_j[v], diag_val[v]);
#endif

    checkStrictDiagonallyDominant(diag_i,diag_j,diag_val, newsize);

    // now quickly generate a test-vector to solve against:

    double *vec = vecallocd(N);
    for(v=0;v<N;v++)
        vec[v]=ran();

    fprintf(stderr,"Left with %d nonzeroes; nonzero density = %lf\n", newsize, newsize/((double)N*N));
    fprintf(stderr,"========== OUTPUTTING ... ==========\n");

    outputSimpleMatrix(newsize, diag_i, diag_j, diag_val, vec);
    fprintf(stderr,"========== THIS IS A DIVIDER ========== \n");
    //outputMatrix(newsize, diag_i, diag_j, diag_val, vec);
    outputMathematicaMatrix(newsize, diag_i, diag_j, diag_val, vec);

    free(vec);
    free(diag_i);
    free(diag_j);
    free(diag_val);

    return 0;
}

/**
 * This function takes some matrix A and produces (almost) A' = A + A^T, which is
 * a symmetric matrix. Note that the diagonals are not added.
 *
 * @return number of nonzeros after transpose has been done.
 */
int addTranspose(int final_size, int* i, int* j, double* v,
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

    int c,c2;
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

int countDiags(int* i, int* j, int nz) {
    int diags=0;
    int c;
    for (c=0; c<nz; c++) {
        if(i[c]==j[c])
            diags++;
    }

    return diags;
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


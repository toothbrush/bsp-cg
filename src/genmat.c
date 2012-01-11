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

    int v;


    for (v=0; v<nz; v++) {

        xs[v]= (double)N * ran();
        ys[v]= (double)N * ran();
        vals[v]= ran()*2-1; // [-1,1]

#ifdef DEBUG
        fprintf(stderr,"generated A[%d][%d]=%lf\n", xs[v],ys[v], vals[v]);
#endif

    }


    // oh my fucking god, look away, it's evil

    int uniques=0;
    bool found = false;

    int* fin_i;
    int* fin_j;
    double* fin_val;

    fin_i = vecalloci(nz);
    fin_j = vecalloci(nz);
    fin_val = vecallocd(nz);
    uniques=0;
    int i;

    // and rebuild an array which actually contains those nonzeroes.
    // this is rather expensive...
    for(v=0;v<nz;v++) {
        // loop through the array, and for each element, check if the preceding 
        // portion of the array contains the same (i,j) pair. count first occurrences
        // of (i,j)'s.

        found=false;
        for(i=0;i<v;i++) {
            if(xs[v] == xs[i] &&
               ys[v] == ys[i]) {
                found = true;
                break;
            }
        }
        if(!found) {
            fin_i[uniques] = xs[v];
            fin_j[uniques] = ys[v];
            fin_val[uniques] = vals[v];

            uniques ++;

        }

    }
    fin_i = realloc(fin_i, uniques*SZINT);
    fin_j = realloc(fin_j, uniques*SZINT);
    fin_val = realloc(fin_val, uniques*SZDBL);

    int diagonals_present = countDiags(fin_i, fin_j, uniques);
#ifdef DEBUG
    fprintf(stderr, "found %d diagonal(s), still need %d more.\n", diagonals_present, (N-diagonals_present));
#endif

    int newsize = uniques + (N - diagonals_present);
    int* diag_i;
    int* diag_j;
    double* diag_val;

    diag_i = vecalloci(newsize);
    diag_j = vecalloci(newsize);
    diag_val = vecallocd(newsize);

    for(v=0; v<uniques; v++) {
        // copy the old values across
        diag_i[v]=fin_i[v];
        diag_j[v]=fin_j[v];
        diag_val[v]=fin_val[v];
    }

    free(fin_i);
    free(fin_j);
    free(fin_val);

    addDiagonal(mu, diag_i, diag_j, diag_val, uniques, diagonals_present, N);
#ifdef DEBUG
    for(v=0;v<newsize;v++)
        fprintf(stderr,"after addDiagonal A[%d][%d]=%lf\n", diag_i[v],diag_j[v], diag_val[v]);

    fprintf(stderr, "Going to make symmetric now...\n");
#endif

    // things must be symmetric, but they aren't, yet
    // ... here's a good place to do the transposing thing.

    int max_symmetric_size = (newsize - N)*2 + N;
    int actual_symmetric_size = -1;
    fin_i = vecalloci(max_symmetric_size);
    fin_j = vecalloci(max_symmetric_size);
    fin_val = vecallocd(max_symmetric_size);

    actual_symmetric_size = addTranspose(newsize,diag_i,diag_j,diag_val,
                              max_symmetric_size,fin_i, fin_j, fin_val);

#ifdef DEBUG
    for(v=0;v<actual_symmetric_size;v++)
        // to make diags stand out.
        if(fin_i[v]==fin_j[v])
            fprintf(stderr,"after transpose A[%d][%d]=%lf \\\\\n", fin_i[v],fin_j[v], fin_val[v]);
        else
            fprintf(stderr,"after transpose A[%d][%d]=%lf\n", fin_i[v],fin_j[v], fin_val[v]);
#endif

    // swap stuff around here:
    free(diag_i); free(diag_j); free(diag_val);

    diag_i = fin_i; diag_j = fin_j; diag_val = fin_val;
    newsize = actual_symmetric_size;

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

    free(xs);
    free(ys);
    free(vals);
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
int addTranspose(int nz, int* i, int* j, double* v,
                 int max,int* out_i, int* out_j, double* out_v) {

    int c,c2;
    int actual_nonzeroes=0;

    int *done_i, *done_j;
    done_i = vecalloci(nz);
    done_j = vecalloci(nz);
    int done_n = 0;

    bool done;

    // also rather expensive
    for(c=0; c<nz; c++) {
        // for each original entry

        if(i[c] == j[c]) {
            // original diagonal

            out_i[actual_nonzeroes] = i[c];
            out_j[actual_nonzeroes] = j[c];
            out_v[actual_nonzeroes] = v[c];
            actual_nonzeroes++;

        } else {
            // original non-diagonal

            // check if it has been done.
            done = false;
            for(c2=0; c2 < done_n; c2++) {
                if((done_i[c2] == i[c] && done_j[c2] == j[c]) ||
                   (done_j[c2] == i[c] && done_i[c2] == j[c])  ) {
                    done = true;
                }
            }

            // if not yet done,
            if(! done) {

                //add v[c] into out_*
                out_i[actual_nonzeroes] = i[c];
                out_j[actual_nonzeroes] = j[c];
                out_v[actual_nonzeroes] = 0;

                out_v[actual_nonzeroes] += v[c];
                for(c2=0; c2 < nz; c2++) {
                    // add all transpose-components

                    if(i[c2] == j[c] &&
                       j[c2] == i[c]) {
                        out_v[actual_nonzeroes] += v[c2];
                    }

                }
                actual_nonzeroes++;

                // and place it in the transposed position.

                out_i[actual_nonzeroes] = j[c];
                out_j[actual_nonzeroes] = i[c];
                out_v[actual_nonzeroes] = out_v[actual_nonzeroes -1];
                actual_nonzeroes++;

                done_i[done_n] = i[c];
                done_j[done_n] = j[c];
                done_n ++;

            }
        }

    }

    free(done_i);
    free(done_j);

    return actual_nonzeroes;

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
        if(i[c] != j[c]) {
            rowtotal[i[c]] += fabs(v[c]);
        }
    }

    // find diagonals:
    for(c = 0; c<nz; c++)
    {
        if(i[c] == j[c]){
            diagonals[i[c]] = v[c];
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

void addDiagonal(double mu, int* i, int* j, double* v, int nz, int diags_present, int diags_needed) {

    int c,c2;
    int pos = nz; // where to start appending.

    bool found;
    // for each diagonal element, check if it exists, if not, append.
    for(c=1;c<=diags_needed;c++)
    {
        found = false;
        for(c2=0;c2<nz;c2++) {
            if(i[c2] == j[c2] && // is a diagonal
               j[c2] == c-1)     // is the diag we are looking for.
            {
                found = true;
#ifdef DEBUG
                fprintf(stderr, "found a diag.\n");
#endif

                v[c2] += mu;
                if (v[c2] <= 0){
                    fprintf(stderr, "ERROR: not all diagonals are >0!\n");
                    exit(22);
                }

                break;

            }
        }
        if (!found) {
            // append
#ifdef DEBUG
            fprintf(stderr, "not found, appending!\n");
#endif
            i[pos] = c-1;
            j[pos] = c-1;
            v[pos] = mu;
            pos++;
        }
    }

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
    printf("somemat // MatrixForm\n");

    printf("\n\n\n");

    fprintf(stderr, "(* ======= vector v follows ====== *)\n");


    printf("vec = {\n");
    // N vector entries, one proc.
    for(c=0;c<N-1;c++) {
        // each line is a value, in order of the vector indices.
        printf("%lf,\n", vec[c]);
    }
    // last line without comma.
    printf("%lf\n};\nvec // MatrixForm\n", vec[N-1]);

    // and finally, for the paranoid:
    printf("\n\nPositiveDefiniteMatrixQ[somemat]\n\nSymmetricMatrixQ[somemat]\n");
    printf("correctAnswer = LinearSolve[somemat,vec];\n");
    printf("correctAnswer // MatrixForm\n");

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


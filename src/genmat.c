#include "libs/vecalloc-seq.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "genmat.h"
#include "libs/paulbool.h"

int N;

int main (int argc, char** argv) {

    enum outputformat out;

    out = MATHEMATICA;

    double sparsity;
    // aim for a nonzero density
    sparsity = 0.02; // nz = sparsity*100% of the size of the matrix

    // read the size of the matrix from command line
    if (argc < 2) {
        printf("Usage: %s N\n", argv[0]);
        exit(-1);
    }

    if(sscanf(argv[1], "%d", &N) != 1) {
        printf("couldn't read command-line argument for N. must be an integer.\n");
        exit(-2);
    }

    double mu;
    mu = 2.0; //scalar for making matrix diagonal-dominant

    int nz = sparsity*N*N;
    int* xs;
    int* ys;
    double* vals;

    fprintf(stderr,"Generating matrix. N=%d, density=%lf, target nz=%d\n",
            N, sparsity, nz);

    //TODO: use random() from stdlib!!!
    srand((unsigned)time(NULL));


    xs = vecalloci(nz);
    ys = vecalloci(nz);
    vals = vecallocd(nz);

    int v;


    for (v=0; v<nz; v++) {

        xs[v]= (double)N * ran();
        ys[v]= (double)N * ran();
        vals[v]= ran()*2-1; // [-1,1]

        fprintf(stderr,"generated A[%d][%d]=%lf\n", xs[v],ys[v], vals[v]);

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
                break;
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
                break;
            }
        }
        if(!found) {
            // copy over.

            fin_i[uniques] = xs[v];
            fin_j[uniques] = ys[v];
            fin_val[uniques] = vals[v];

            uniques ++;

        }

    }

    /*
    addTranspose(fin_i, fin_j, fin_val, uniques);


    for(v=0;v<uniques;v++)
        // to make diags stand out.
        if(fin_i[v]==fin_j[v])
            fprintf(stderr,"after transpose A[%d][%d]=%lf \\\\\n", fin_i[v],fin_j[v], fin_val[v]);
        else
            fprintf(stderr,"after transpose A[%d][%d]=%lf\n", fin_i[v],fin_j[v], fin_val[v]);

     */

    int diagonals_present = countDiags(fin_i, fin_j, uniques);
    fprintf(stderr, "found %d diagonal(s), still need %d more.\n", diagonals_present, (N-diagonals_present));

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

    addDiagonal(mu, diag_i, diag_j, diag_val, uniques, diagonals_present, N);
    for(v=0;v<newsize;v++)
        fprintf(stderr,"after addDiagonal A[%d][%d]=%lf\n", diag_i[v],diag_j[v], diag_val[v]);

    // things must be symmetric, but they aren't FIXME
    // ... here's a good place to do the transposing thing.

    checkStrictDiagonallyDominant(diag_i,diag_j,diag_val, newsize);

    fprintf(stderr,"Left with %d nonzeroes; nonzero density = %lf\n", uniques, newsize/((double)N*N));
    fprintf(stderr,"========== OUTPUTTING ... ==========\n");

    if(out == SIMPLE) {
        outputSimpleMatrix(newsize, diag_i, diag_j, diag_val);
    } else if (out == EMM) {
        outputMatrix(newsize, diag_i, diag_j, diag_val);
    } else if (out == MATHEMATICA) {
        outputMathematicaMatrix(newsize, diag_i, diag_j, diag_val);
    }

    free(xs);
    free(ys);
    free(vals);
    free(fin_i);
    free(fin_j);
    free(fin_val);
    free(diag_i);
    free(diag_j);
    free(diag_val);

    return 0;
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
                fprintf(stderr, "found a diag.\n");

                v[c2] += mu;
                if (v[c2] <= 0)
                    fprintf(stderr, "ERROR: not all diagonals are >0!\n");

                break;

            }
        }
        if (!found) {
            // append
            fprintf(stderr, "not found, appending!\n");
            i[pos] = c-1;
            j[pos] = c-1;
            v[pos] = mu;
            pos++;
        }
    }

}

// FIXME: this is still dodgy:
void addTranspose(int* i, int* j, double* v, int nz) {


    int c,c2;
    int tx;

    // we have to cache already-added (i,j)'s

    int twiddled=0;
    int* done_i;
    int* done_j;
    done_i = vecalloci(nz);
    done_j = vecalloci(nz);

    bool already_done;
    for(c=0; c<nz; c++) {
        // find transpose, if it exists, add it to current.

        for(tx=0; tx<nz; tx++) {

            if (i[tx] == j[c] &&
                j[tx] == i[c] &&
                c     != tx      // don't double everything!
                ) {

                already_done = false;
                for(c2=0; c2<twiddled; c2++) {
                    if (i[tx] == done_i[c2] &&
                        j[tx] == done_j[c2]) {
                        already_done = true;
                        break;
                    }
                }

                if(already_done) {
                    v[c] = v[tx]; // just copy, the add has already been done.

                }
                else
                {
                    v[c] += v[tx];
                    done_i[twiddled] = i[c];
                    done_j[twiddled] = j[c];
                    twiddled++;
                    fprintf(stderr,"transpose of (%d,%d)!\n", i[tx], j[tx]);
                }

            }


        }

    }
    free(done_i);
    free(done_j);
    fprintf(stderr,"twiddled = %d\n", twiddled);
}


double ran() {

    return ((double)rand()/(double)RAND_MAX);

}

void outputMathematicaMatrix(int nz, int*i, int*j, double*v) {

    // here we'll print the matrix in Mathematica format

    // Mathematica "header"

    printf("somemat = SparseArray[ { \n");
    // the value lines: i j value:
    int c;
    for(c=0;c<nz-1;c++) {

        printf("{%d, %d} -> %lf,\n", i[c]+1, j[c]+1, v[c]); // Mondriaan expects 1-based coordinates.

    }
    printf("{%d, %d} -> %lf\n", i[nz-1]+1, j[nz-1]+1, v[nz-1]); // Mondriaan expects 1-based coordinates.

    printf("} ] ;\n");
    printf("somemat // MatrixForm\n");

    printf("\n\n\n");

    fprintf(stderr, "======= vector v follows ======\n");


    printf("vec = {\n");
    // N vector entries, one proc.
    for(c=0;c<N-1;c++) {
        // each line is a value, in order of the vector indices.
        //printf("%d %d %lf\n", c+1, 1, ran());
        printf("%lf,\n", ran());
    }
    // last line without comma.
    printf("%lf\n};\nvec // MatrixForm\n", ran());

    // and finally, for the paranoid:
    printf("\n\nPositiveDefiniteMatrixQ[somemat]\n\nSymmetricMatrixQ[somemat]\n");
    printf("correctAnswer = LinearSolve[somemat,vec]\n");

}
void outputSimpleMatrix(int nz, int*i, int*j, double*v) {

    // here we'll print the matrix in simple format, for one proc.

    // no header.
    //size line: m n nz
    printf("%d %d %d %d\n", N, N, nz, 1); // one processor, i.e. not distributed
    // begin and end bounds of proc 1 are 0 and nz
    printf("0\n%d\n", nz);

    // next the value lines: i j value:
    int c;
    for(c=0;c<nz;c++) {

        printf("%d %d %lf\n", i[c]+1, j[c]+1, v[c]); // Mondriaan expects 1-based coordinates.

    }

    fprintf(stderr, "======= vector v follows ======\n");


    // N vector entries, one proc.
    printf("%d %d\n", N, 1);
    for(c=0;c<N;c++) {
        // each line is
        //   coordinate processor value
        printf("%d %d %lf\n", c+1, 1, ran());
    }


}
void outputMatrix(int nz, int*i, int*j, double*v) {

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
        printf("%lf\n", ran());
    }


}


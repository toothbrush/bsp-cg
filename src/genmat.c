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

            uniques ++;

        }

    }

    addTranspose(fin_i, fin_j, fin_val, uniques);


    for(v=0;v<uniques;v++)
        printf("left with A[%d][%d]=%lf\n", fin_i[v],fin_j[v], fin_val[v]);



    printf("Generated %d nonzeroes\n", uniques);



    free(xs);
    free(ys);
    free(vals);
    free(fin_i);
    free(fin_j);
    free(fin_val);

    return 0;
}

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
                c != tx // don't double everything!
                ) {

                already_done = false;
                for(c2=0;c2<twiddled;c2++) {
                    if (i[tx] == done_i[c2] &&
                        j[tx] == done_j[c2]) {
                        already_done = true;
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
                    printf("transpose of (%d,%d)!\n", i[tx], j[tx]);
                }

            }


        }

    }
    free(done_i);
    free(done_j);
    printf("twiddled = %d\n", twiddled);
}


double ran() {

    return ((double)rand()/(double)RAND_MAX);

}

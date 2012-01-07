void   checkStrictDiagonallyDominant(int* i, int* j, double* v, int nz);
int    countDiags(int* i, int* j, int nz);
void   addDiagonal(double mu, int* i, int* j, double* v, int nz, int diags_present, int diags_needed);
int    addTranspose(int nz, int* i, int* j, double* v, int* out_i, int* out_j, double* out_v);
double ran();
void   outputMatrix(int nz, int*i, int*j, double*v);
void   outputMathematicaMatrix(int nz, int*i, int*j, double*v);
void   outputSimpleMatrix(int nz, int*i, int*j, double*v);


enum outputformat {
    EMM,
    SIMPLE,
    MATHEMATICA
};



#include "libs/paulbool.h"

void   checkStrictDiagonallyDominant(int* i, int* j, double* v, int nz);
int    countDiags(int* i, int* j, int nz);
void   addTranspose(int nz, int* i, int* j, double* v, int maxsize);
double ran();
void addDiagonal(double mu, int* i, int* j, double* v, int nz, int diags_needed, bool* diags_done);
// output functions
void   outputMatrix(int nz, int*i, int*j, double*v, double*vec);
void   outputMathematicaMatrix(int nz, int*i, int*j, double*v, double*vec);
void   outputSimpleMatrix(int nz, int*i, int*j, double*v, double*vec);

#include "libs/paulbool.h"

void   checkStrictDiagonallyDominant(int* i, int* j, double* v, int nz);
int    countDiags(int* i, int* j, int nz);
int    addTranspose(int nz, int* i, int* j, double* v, int maxsize);
double ran();

// output functions
void   outputMatrix(int nz, int*i, int*j, double*v, double*vec);
void   outputMathematicaMatrix(int nz, int*i, int*j, double*v, double*vec);
void   outputSimpleMatrix(int nz, int*i, int*j, double*v, double*vec);

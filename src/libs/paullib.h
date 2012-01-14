#include "paulbool.h"

bool file_exists(const char * filename);
void out(int proc, char*at, const char *fmt, ...);
void negate(int n, double* v);
void scalevec(int n, double factor, double*vec);
void local_axpy(int n, double a, double* x, double* y,double* result);
void zero(int nv, double* v);
void one(int nv, double* v);

double ran();

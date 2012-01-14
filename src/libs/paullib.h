#include "paulbool.h"

bool file_exists(const char * filename);
void out(int proc, char*at, const char *fmt, ...);
void negate(int n, double* v);
void copyvec(int nv, int nu, double* v, double* u, int* procu, int* indu);
void scalevec(int n, double factor, double*vec);
void addvec(int n, double* dest, double* a, double*b);
void local_axpy(int n, double a, double* x, double* y,double* result);
void zero(int nv, double* v);
void one(int nv, double* v);

double ran();

#include "paulbool.h"

bool file_exists(const char * filename);
void out(int proc, char*at, const char *fmt, ...);
void negate(int n, double* v);
void copyvec(int n, double* dest, double* src);
void scalevec(int n, double factor, double*vec);
void addvec(int n, double* dest, double* a, double*b);
void axpy(int nv, double a, double* x, double* y,double* result);
void zero(int nv, double* v);
void one(int nv, double* v);

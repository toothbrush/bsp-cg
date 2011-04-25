/*
 * Author: Paul van der Walt
 *
 * March 2011
 *
 * Library with a few useful functions.
 */

#include <stdio.h>
#include <stdarg.h>

void negate(int n, double* v)
{
    int i;
    for(i = 0; i<n; i++)
        v[i] *= -1;
    
}

void copyvec(int n, double* dest, double* src)
{
    int i;
    for(i = 0; i<n; i++)
        dest[i] = src[i];
    
}
void scalevec(int n, double factor, double*vec)
{
    int i;
    for(i = 0; i<n; i++)
        vec[i] *= factor;
}
void addvec(int n, double* dest, double* a, double*b)
{
    int i;
    for(i = 0; i<n; i++)
        dest[i] = a[i] + b[i];
}

/*
 * This function wraps printf to redirect debug
 * output to stderr. This is useful if we (later) want
 * to create a tool pipeline like
 *
 * $ gen-matrix | proc-matrix | print-matrix
 *
 * or whatever. 
 */
void out(int proc, char*at, const char *fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);

    char extended_fmt[1024];
    sprintf(extended_fmt, "%s (P%d) -- \t%s", at, proc, fmt);
    vfprintf(stderr, extended_fmt, argp);

    va_end(argp);
}

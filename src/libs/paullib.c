/*
 * Author: Paul van der Walt
 *
 * March 2011
 *
 * Library with a few useful functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "paullib.h"

/*
 * useful function to stat a file
 */
bool file_exists(const char * filename)
{
    FILE *file;
    if((file = fopen(filename, "r")))
    {
        fclose(file);
        return true;
    }
    return false;
}

/*
 * negate all of a vector's components
 */
void negate(int n, double* v)
{
    int i;
    for(i = 0; i<n; i++)
        v[i] *= -1.0;

}

/*
 * copy vec src into dest
 */
void copyvec(int n, double* src,
                    double* dest)
{
    int i;
    for(i = 0; i<n; i++) {
        dest[i] = src[i];
    }

}

/*
 * make all vector components 0
 */
void zero (int nv, double * a)
{
    int i;
    for (i = 0; i < nv; i++)
        a[i] = 0.0;
}

/*
 * Calculate ax+y.
 */
void axpy (int nv, double a, double* x, double* y,double* result) {

    int i;
    for (i = 0; i< nv; i++) {
        result[i] = a * x[i] + y[i];
    }

}

/*
 * Scale a vector vec by a factor.
 */
void scalevec(int n, double factor, double*vec)
{
    int i;
    for(i = 0; i<n; i++)
        vec[i] *= factor;
}

/*
 * Add a and b and place result in dest. Assumes dest
 * is malloc'd already.
 */
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
    // may be a bit hacky, but this gives different colours to different proc's output.
    // only possible on POSIX I'm afraid. Untested on anything but Linux+zsh as yet.
    // Turns out OS X also does the trick.
    if(proc == -1) {
        sprintf(extended_fmt, "%s (P.) -- \t%s", at, fmt);
    } else {
        // on Huygens, using mpcc, we can't use colours, so skip them if we're
        // not using GCC.
#ifdef __GNUC__
        sprintf(extended_fmt, "\e[0;%dm%s (P%d) -- \t%s\e[0m", 31+proc, at, proc, fmt);
#else
        sprintf(extended_fmt, "%s (P%d) -- \t%s", at, proc, fmt);
#endif
    }

#ifdef DEBUG
    vfprintf(stderr, extended_fmt, argp);
    fflush(stderr);
#endif

    va_end(argp);
}

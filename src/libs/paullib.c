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
#include "bspedupack.h"

/*
 * return a random double in the interval [0,1)
 */
double ran() {

    return ((double)random()/(double)RAND_MAX);

}


/*
 * useful function to stat a file for existence
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
 * make all vector components 0
 */
void zero (int nv, double * a)
{
    int i;
    for (i = 0; i < nv; i++)
        a[i] = 0.0;
}

/*
 * Calculate ax+y for all local vector components, and
 * store the result in array "result".
 */
void local_axpy (int n, double a, double* x, double* y,double* result) {

    int i;
    for (i = 0; i< n; i++) {
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
 * This function wraps printf to redirect debug
 * output to stderr. This is useful if we want
 * to filter output, like for example
 *
 * $ ./bin/cg 2> /dev/null
 *
 * or whatever.
 */
void out(int proc, char*at, const char *fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);

    char extended_fmt[1024];
    // may be a bit hacky, but this gives different colours to different proc's
    // output.  only possible on POSIX I'm afraid. Untested on anything but
    // Linux+{zsh,bash} and OS X.
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

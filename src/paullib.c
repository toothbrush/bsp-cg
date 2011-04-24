/*
 * Author: Paul van der Walt
 *
 * March 2011
 *
 * Library with a few useful functions.
 */

#include <stdio.h>
#include <stdarg.h>

void negate(double* v)
{

}

void copyvec(double* dest, double* src)
{

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
void out(const char *fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);

    vfprintf(stderr, fmt, argp);

    va_end(argp);
}

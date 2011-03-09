#include <stdio.h>
#include <stdarg.h>

extern char *itoa(int, char *, int);

void out(const char *fmt, ...)
{
    va_list argp;

    va_start(argp, fmt);

    printf("%s", fmt);
    printf(fmt, argp);

    va_end(argp);
}

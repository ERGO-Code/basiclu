/*
 * lu_timer.c
 *
 * Wall clock timer copied from T. Davis, SuiteSparse.
 *
 */

#ifndef BASICLU_NOTIMER
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#endif
#include "lu_timer.h"

void lu_tic (double tic[2])
{
#ifndef BASICLU_NOTIMER
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    tic[0] = (double) t.tv_sec;
    tic[1] = (double) t.tv_nsec;
#else
    tic[0] = 0.0;
    tic[1] = 0.0;
#endif
}

double lu_toc (const double tic[2])
{
    double toc[2];
    lu_tic(toc);
    return (toc[0] - tic[0]) + 1e-9*(toc[1] - tic[1]);
}

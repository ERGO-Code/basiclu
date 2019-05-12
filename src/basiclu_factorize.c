/*
 * basiclu_factorize.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "lu_internal.h"

lu_int basiclu_factorize
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
    const lu_int Bbegin[],
    const lu_int Bend[],
    const lu_int Bi[],
    const double Bx[],
    lu_int c0ntinue
)
{
    return basiclu_factorize_buckets(istore, xstore, Li, Lx, Ui, Ux, Wi, Wx,
                                     Bbegin, Bend, Bi, Bx, NULL, c0ntinue);
}

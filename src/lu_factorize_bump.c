/* 
 * lu_factorize_bump.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Bump factorization driver routine
 *
 */

#include "lu_internal.h"

lu_int lu_factorize_bump(struct lu *this)
{
    const lu_int m  = this->m;
    lu_int *pinv    = this->pinv;
    lu_int *qinv    = this->qinv;

    lu_int pivot_row, pivot_col;
    lu_int status = BASICLU_OK;

    while (this->rank + this->rankdef < m)
    {
        /*
         * Find pivot element. Markowitz search need not be called if the
         * previous call to lu_pivot() returned for reallocation. In this case
         * this->pivot_col is valid.
         */
        if (this->pivot_col < 0)
            lu_markowitz(this);
        assert(this->pivot_col >= 0);
        pivot_row = this->pivot_row;
        pivot_col = this->pivot_col;
        assert(pinv[pivot_row] == -1);
        assert(qinv[pivot_col] == -1);

        /* Eliminate pivot. This may require reallocation. */
        status = lu_pivot(this);
        if (status != BASICLU_OK)
            break;

        pinv[pivot_row] = this->rank;
        qinv[pivot_col] = this->rank;
        this->pivot_col = -1;
        this->pivot_row = -1;
        this->rank++;
    }

    return status;
}

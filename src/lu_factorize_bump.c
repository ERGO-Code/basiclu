/* 
 * lu_factorize_bump.c
 *
 * Copyright (C) 2016-2019  ERGO-Code
 *
 * Bump factorization driver routine
 *
 */

#include "lu_internal.h"
#include "lu_list.h"

static void lu_next_bucket(struct lu *this);

lu_int lu_factorize_bump(struct lu *this)
{
    const lu_int m          = this->m;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *pinv            = this->pinv;
    lu_int *qinv            = this->qinv;
    lu_int status = BASICLU_OK;

    while (this->rank + this->rankdef < m)
    {
        if (this->bucket_ptr == this->rank+this->rankdef+1)
            lu_next_bucket(this);
        assert(this->ncol_active > 0);

        /*
         * Find pivot element. Markowitz search need not be called if the
         * previous call to lu_pivot() returned for reallocation. In this case
         * this->pivot_col is valid.
         */
        if (this->pivot_col < 0)
            lu_markowitz(this);
        assert(this->pivot_col >= 0);

        if (this->pivot_row < 0)
        {
            /* Eliminate empty column without choosing a pivot. */
            lu_list_remove(colcount_flink, colcount_blink, this->pivot_col);
            this->pivot_col = -1;
            this->rankdef++;
            this->ncol_active--;
        }
        else
        {
            /* Eliminate pivot. This may require reallocation. */
            assert(pinv[this->pivot_row] == -1);
            assert(qinv[this->pivot_col] == -1);
            status = lu_pivot(this);
            if (status != BASICLU_OK)
                break;
            pinv[this->pivot_row] = this->rank;
            qinv[this->pivot_col] = this->rank;
            this->pivot_col = -1;
            this->pivot_row = -1;
            this->rank++;
            this->ncol_active--;
        }
    }
    return status;
}

static void lu_next_bucket(struct lu *this)
{
    const lu_int m          = this->m;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    const lu_int *Wbegin    = this->Wbegin;
    const lu_int *Wend      = this->Wend;
    const lu_int *Lbegin_p  = this->Lbegin_p;
    lu_int bucket_ptr       = this->bucket_ptr;
    lu_int j, nz;

    assert(bucket_ptr >= 1 && bucket_ptr <= m);
    assert(Lbegin_p[bucket_ptr] < 0);
    assert(this->ncol_active == 0);

    do {
        j = Lbegin_p[bucket_ptr++];
        if (j < 0) j = -j-1;    /* must unflip first index */
        assert(this->qinv[j] < 0);
        nz = Wend[j] - Wbegin[j];
        lu_list_add(j, nz, colcount_flink, colcount_blink, m, &this->min_colnz);
        this->ncol_active++;
    } while (bucket_ptr <= m && Lbegin_p[bucket_ptr] >= 0);

    this->bucket_ptr = bucket_ptr;
}

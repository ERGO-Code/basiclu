/*
 * basiclu_obj_maxvolume.c
 *
 * Copyright (C) 2016-2017  ERGO-Code
 *
 */

#include "lu_internal.h"

/*
 * factorize() - factorize A[:,basis]
 */
static lu_int factorize(struct basiclu_object *obj,
                        const lu_int *Ap,
                        const lu_int *Ai,
                        const double *Ax,
                        const lu_int *basis)
{
    double *xstore = obj->xstore;
    const lu_int m = xstore[BASICLU_DIM];
    lu_int *begin = NULL;
    lu_int *end = NULL;
    lu_int i, status = BASICLU_OK;

    begin = malloc(m*sizeof(lu_int));
    end = malloc(m*sizeof(lu_int));
    if (!begin || !end) {
        status = BASICLU_ERROR_out_of_memory;
        goto cleanup;
    }
    for (i = 0; i < m; i++)
    {
        begin[i] = Ap[basis[i]];
        end[i] = Ap[basis[i]+1];
    }

    status = basiclu_obj_factorize(obj, begin, end, Ai, Ax);

cleanup:
    if (begin) free(begin);
    if (end) free(end);

    return status;
}

/*
 * refactorize_if_needed() - refactorize the basis if required or favourable
 *
 * The basis matrix is refactorized if
 * - the maximum number of updates is reached, or
 * - the previous update had a large pivot error, or
 * - it is favourable for performance
 *
 * factorize() is called for the actual factorization.
 *
 * Note: refactorize_if_needed() will not do an initial factorization.
 */
static lu_int refactorize_if_needed(struct basiclu_object *obj,
                                    const lu_int *Ap,
                                    const lu_int *Ai,
                                    const double *Ax,
                                    const lu_int *basis)
{
    lu_int status = BASICLU_OK;
    const double piverr_tol = 1e-8;
    double *xstore = obj->xstore;

    if (xstore[BASICLU_NFORREST] == xstore[BASICLU_DIM] ||
        xstore[BASICLU_PIVOT_ERROR] > piverr_tol ||
        xstore[BASICLU_UPDATE_COST] > 1.0)
        status = factorize(obj, Ap, Ai, Ax, basis);
    return status;
}

/*
 * basiclu_obj_maxvolume() - one pass over columns of A doing basis updates
 *
 * For each column a_j not in B, compute lhs = B^{-1}*a_j and find the maximum
 * entry lhs[imax]. If it is bigger than @volumetol in absolute value, then
 * replace position imax of the basis by index j. On return *p_nupdate is the
 * number of basis updates done.
 */
lu_int basiclu_obj_maxvolume
(
    struct basiclu_object *obj,
    lu_int ncol,
    const lu_int Ap[],
    const lu_int Ai[],
    const double Ax[],
    lu_int basis[],
    lu_int isbasic[],
    double volumetol,
    lu_int *p_nupdate
)
{
    lu_int i, j, k;
    lu_int nzrhs, imax, begin, nupdate = 0;
    double xtbl, xmax;
    lu_int status = BASICLU_OK;

    if (volumetol < 1.0)
    {
        status = BASICLU_ERROR_invalid_argument;
        goto cleanup;
    }

    /* Compute initial factorization. */
    status = factorize(obj, Ap, Ai, Ax, basis);
    if (status != BASICLU_OK)
        goto cleanup;

    for (j = 0; j < ncol; j++)
    {
        if (isbasic[j])
            continue;

        /* compute B^{-1}*a_j */
        nzrhs = Ap[j+1] - Ap[j];
        begin = Ap[j];
        status = basiclu_obj_solve_for_update(obj, nzrhs, Ai+begin, Ax+begin,
                                              'N', 1);
        if (status != BASICLU_OK)
            goto cleanup;

        /* Find the maximum entry. */
        xmax = 0.0;
        xtbl = 0.0;
        imax = 0;
        for (k = 0; k < obj->nzlhs; k++) {
            i = obj->ilhs[k];
            if (fabs(obj->lhs[i]) > xmax) {
                xtbl = obj->lhs[i];
                xmax = fabs(xtbl);
                imax = i;
            }
        }

        if (xmax <= volumetol)
            continue;

        /* Update basis. */
        isbasic[basis[imax]] = 0;
        isbasic[j] = 1;
        basis[imax] = j;
        nupdate++;

        /* Prepare to update factorization. */
        status = basiclu_obj_solve_for_update(obj, 0, &imax, NULL, 'T', 0);
        if (status != BASICLU_OK)
            goto cleanup;

        status = basiclu_obj_update(obj, xtbl);
        if (status != BASICLU_OK)
            goto cleanup;

        status = refactorize_if_needed(obj, Ap, Ai, Ax, basis);
        if (status != BASICLU_OK)
            goto cleanup;
    }

cleanup:
    if (p_nupdate) *p_nupdate = nupdate;
    return status;
}

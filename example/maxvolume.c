/*
 * maxvolume.c
 *
 * Copyright (C) 2016-2017  Lukas Schork
 *
 * Command line program to compute a maximum volume basis [1,2]
 *
 * Purpose:
 *
 * Given an m-by-n matrix A, usually m < n, find m columns of A that form a
 * nonsingular matrix B with maximum volume in A. When A does not have m
 * linearly independent columns (in particular when m > n), then the remaining
 * columns of B are columns of the identity matrix. The number of such columns
 * is the row rank deficiency of A.
 *
 * The basis B has (local) maximum volume in A iff for every column a_j of A
 * all entries in B^{-1}*a_j are less than or equal to 1 in absolute value.
 *
 * Algorithm:
 *
 * Append 1e-8 * identity matrix to A and use these columns as initial basis.
 * Then repeatedly pass over the columns of A. If column a_j is not in B and
 * B^{-1}*a_j has an entry larger than @volumetol in absolute value, then a_j
 * replaces one column in B. The algorithm stops when one pass over the columns
 * of A did not change B, or when @maxpass passes are done.
 *
 * Usually the major cost of the algorithm is computing B^{-1}*a_j. When these
 * vectors are dense, then one pass over the columns of A has time complexity
 * Omega(m*n) and the algorithm is impractical for large scale matrices. When
 * these vectors are sparse (as they frequently are in linear programming
 * problems), then the BASICLU routines will exploit the sparsity and the
 * computational cost for one pass is proportional to the number of nonzeros in
 * B^{-1}*A. See the LP problems in data/.
 *
 * Note:
 *
 * On Maragal_4 the algorithm fails with default parameters because a
 * refactorization reports a singularity in the basis. Debugging has shown
 * that the final entry in the active submatrix is 2.03e-15, which is below
 * the default absolute pivot tolerance (1e-14). Changing the absolute pivot
 * tolerance to 1e-15 makes the algorithm complete, as does tightening the
 * relative pivot tolerance to 0.5. In the latter case I'm not sure if that
 * results accidently from a different refactorization point or if it fixes
 * numerical instability in the factorization.
 *
 * Parameters:
 *
 * maxvolume <matrix.mtx> [volumetol [maxpass]]
 *
 * <matrix.mtx>: name of a file containing the matrix A in matrix market format
 *               (real, sparse, general)
 * volumetol:    tolerance >= 1 on the absolute value of B^{-1}*a_j
 *               default: 1.1
 * maxpass:      maximum # passes over the columns of A
 *               default: 2
 *
 * [1] C. T. Pan, "On the existence and computation of rank-revealing LU
 *     factorizations". Linear Algebra Appl., 316(1-3), pp. 199-222, 2000
 *
 * [2] S. A. Goreinov, I. V. Oseledets, D. V. Savostyanov, E. E. Tyrtyshnikov,
 *     N. L. Zamarashkin, "How to find a good submatrix". In "Matrix methods:
 *     theory, algorithms and applications", pp. 247-256. World Sci. Publ.,
 *     Hackensack, NJ, 2010.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "basiclu.h"
#include "mmio.h"

/* internal error codes, supplementing BASICLU status codes */
#define IO_ERROR 102

int main(int argc, const char *argv[])
{
    lu_int *Ap      = NULL;
    lu_int *Ai      = NULL;
    double *Ax      = NULL;
    lu_int *basis   = NULL;
    lu_int *isbasic = NULL;
    lu_int *count   = NULL;
    int *mm_I       = NULL;
    int *mm_J       = NULL;
    double *mm_val  = NULL;
    int mm_M, mm_N, mm_nz;
    lu_int m, n, nz, i, j, l, put, err = 0;
    struct basiclu_object factor;

    double volumetol = 1.1;     /* command line parameters */
    long maxpass = 2;

    basiclu_obj_initialize(&factor, 0); /* nullify members */

    if (argc < 2 || argc > 4) {
        printf(" usage: maxvolume <matrix.mtx> [volumetol [maxpass]]\n");
        return 101;
    }
    if (argc >= 3)
        volumetol = atof(argv[2]);
    if (argc >= 4)
        maxpass = atol(argv[3]);
    volumetol = fmax(volumetol, 1.0);

    err = mm_read_unsymmetric_sparse(argv[1], &mm_M, &mm_N, &mm_nz,
                                     &mm_val, &mm_I, &mm_J);
    if (err) {
        err = IO_ERROR;
        goto cleanup;
    }
    printf(" matrix: %d rows, %d columns, %d nonzeros\n", mm_M, mm_N, mm_nz);
    printf(" parameters: volumetol = %.2f, maxpass = %ld\n",
           volumetol, maxpass);

    m = mm_M;                   /* convert to lu_int */
    n = mm_N;
    nz = mm_nz;

    /*
     * Convert coordinate format to compressed column format and append 1e-8 *
     * identity matrix. Use work array @count to count the number of nonzeros
     * per column. While filling the matrix, @count[j] is the next unused
     * position in column j.
     */
    Ap = malloc((m+n+1)*sizeof(lu_int));
    Ai = malloc((nz+m)*sizeof(lu_int));
    Ax = malloc((nz+m)*sizeof(double));
    count = calloc(n, sizeof(lu_int));
    if (!Ap || !Ai || !Ax || !count) {
        err = BASICLU_ERROR_out_of_memory;
        goto cleanup;
    }

    for (l = 0; l < nz; l++)
        count[mm_J[l]]++;
    put = 0;
    for (j = 0; j < n; j++)
    {
        Ap[j] = put;
        put += count[j];
        count[j] = Ap[j];
    }
    Ap[n] = put;
    for (l = 0; l < nz; l++)
    {
        j = mm_J[l];
        put = count[j]++;
        Ai[put] = mm_I[l];
        Ax[put] = mm_val[l];
    }
    for (i = 0; i < m; i++)
    {
        Ai[Ap[n+i]] = i;
        Ax[Ap[n+i]] = 1e-8;
        Ap[n+i+1] = Ap[n+i] + 1;
    }
    free(mm_I);
    free(mm_J);
    free(mm_val);
    free(count);
    mm_I = NULL;
    mm_J = NULL;
    mm_val = NULL;
    count = NULL;

    /*
     * Initialize @factor. Initialize @basis, @isbasic to logical basis and run
     * maxvolume() until the basis does not change any more or @maxpass passes
     * through the matrix are done.
     */
    err = basiclu_obj_initialize(&factor, m);
    if (err != BASICLU_OK)
        goto cleanup;

    basis = malloc(m*sizeof(lu_int));
    isbasic = calloc(m+n, sizeof(lu_int));
    if (!basis || !isbasic) {
        err = BASICLU_ERROR_out_of_memory;
        goto cleanup;
    }
    for (i = 0; i < m; i++)
        basis[i] = n+i;
    for (j = n; j < m+n; j++)
        isbasic[j] = 1;

    long pass, changed = 1;
    for (pass = 0; pass < maxpass && changed; pass++)
    {
        lu_int nupdate;
        err = basiclu_obj_maxvolume(&factor, m+n, Ap, Ai, Ax, basis, isbasic,
                                    volumetol, &nupdate);
        if (err != BASICLU_OK)
            goto cleanup;
        changed = nupdate > 0;
        printf(" pass %d: %d updates\n", pass+1, nupdate);
    }

    /*
     * The number of logical columns in the basis is the number of row rank
     * deficiencies in the input matrix. Print total statistics.
     */
    long rankdef = 0;
    for (j = n; j < m+n; j++)
        rankdef += isbasic[j] != 0;

    long nupdate = factor.xstore[BASICLU_NUPDATE_TOTAL];
    long nforrest = factor.xstore[BASICLU_NFORREST_TOTAL];
    long nperm = nupdate-nforrest;
    long nfactorize = factor.xstore[BASICLU_NFACTORIZE];
    double time_factorize = factor.xstore[BASICLU_TIME_FACTORIZE_TOTAL];
    double time_solve = factor.xstore[BASICLU_TIME_SOLVE_TOTAL];
    double time_update = factor.xstore[BASICLU_TIME_UPDATE_TOTAL];

    printf(" status:               %s\n",
           changed ? "max # passes done" : "optimal basis");
    printf(" # passes:             %ld\n", pass);
    printf(" row rank deficiency:  %ld\n", rankdef);
    printf(" updates [perm + FT]:  %ld [%ld + %ld]\n",
           nupdate, nperm, nforrest);
    printf(" # factorizations:     %ld\n", nfactorize);
    printf(" time factorize:       %.2f sec\n", time_factorize);
    printf(" time solve:           %.2f sec\n", time_solve);
    printf(" time update:          %.2f sec\n", time_update);

cleanup:
    if (Ap) free(Ap);
    if (Ai) free(Ai);
    if (Ax) free(Ax);
    if (basis) free(basis);
    if (isbasic) free(isbasic);
    if (count) free(count);
    if (mm_I) free(mm_I);
    if (mm_J) free(mm_J);
    if (mm_val) free(mm_val);
    basiclu_obj_free(&factor);

    if (err != BASICLU_OK)
        printf(" error (%ld)\n", (long) err);
    return err;
}

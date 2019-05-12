lu_int basiclu_factorize_buckets
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
    const lu_int *buckets,
    lu_int c0ntinue
);

/*
Purpose:

    Performs the same task as basiclu_factorize(), just that the column ordering
    of the factorization is restricted to pivot on all columns in bucket 0,
    before those in bucket 1, and so on. If buckets == NULL, then no restriction
    is imposed and the function performs identically to basiclu_factorize().

Return:

    see basiclu_factorize().

Arguments:

    const lu_int *buckets

        buckets[j] holds the bucket to which column j belongs. Buckets are
        numbered from 0 to m-1. All columns in bucket 0 are chosen as pivot
        column before those in bucket 1, and so on. When a column does not
        contain an eligible pivot element (by means of the absolute pivot
        tolerance), then it is ordered last in the pivot sequence, outside its
        bucket. When buckets == NULL, then the function performs identically to
        basiclu_factorize().

    for the remaining arguments see basiclu_factorize().

Parameters:

    xstore[BASICLU_SEARCH_ROWS]

        This parameter has no effect when buckets are given, in which case only
        columns are searched for pivot elements.

    for the remaining parameters see basiclu_factorize().

Info:

    xstore[BASICLU_STATUS]: status code.

        BASICLU_ERROR_invalid_argument

            The matrix is invalid (a column has a negative number of entries,
            a row index is out of range, or a column has duplicate entries) or
            buckets is invalid (buckets[j] < 0 or buckets[j] >= m for some j).

        for the remaining status codes see basiclu_factorize().

    xstore[BASICLU_NBUCKETS] number of (nonempty) buckets. 0 if buckets == NULL.

    for the remaining info fields see basiclu_factorize().
*/

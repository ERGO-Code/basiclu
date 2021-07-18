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
);

/*
Purpose:

    Make one pass over the columns of a rectangular (ncol >= nrow) matrix and
    pivot each nonbasic column into the basis when it increases the volume (i.e.
    the absolute value of the determinant) of the basis matrix. This is one main
    loop of the "maximum volume" algorithm described in [1,2].

    [1] C. T. Pan, "On the existence and computation of rank-revealing LU
        factorizations". Linear Algebra Appl., 316(1-3), pp. 199-222, 2000

    [2] S. A. Goreinov, I. V. Oseledets, D. V. Savostyanov, E. E. Tyrtyshnikov,
        N. L. Zamarashkin, "How to find a good submatrix". In "Matrix methods:
        theory, algorithms and applications", pp. 247-256. World Sci. Publ.,
        Hackensack, NJ, 2010.

Return:

    BASICLU_ERROR_invalid_argument when volumetol is less than 1.0.
    BASICLU_ERROR_out_of_memory when memory allocation in this function failed.

    The return code from a basiclu_obj_* function called when not BASICLU_OK.
    (Note that BASICLU_WARNING_singular_matrix means that the algorithm failed.)

    BASICLU_OK otherwise.

Arguments:

    struct basiclu_object *obj

        Pointer to an initialized BASICLU object. The dimension of the object
        specifies the number of rows of the matrix.

    lu_int ncol
    const lu_int Ap[ncol+1]
    const lu_int Ai[]
    const double Ax[]

        Matrix A in compressed sparse column format. Column j contains elements

            Ai[Ap[j] .. Ap[j+1]-1], Ax[Ap[j] .. Ap[j+1]-1].

        The columns must not contain duplicate row indices. The row indices per
        column need not be sorted.

   lu_int basis[nrow]

        On entry holds the column indices of A that form the initial basis. On
        return holds the updated basis. A basis defines a square nonsingular
        submatrix of A. If the initial basis is (numerically) singular, then the
        initial LU factorization will fail and BASICLU_WARNING_singular_matrix
        is returned.

   lu_int isbasic[ncol]

        This array must be consistent with basis[] on entry, and is consistent
        on return. isbasic[j] must be nonzero iff column j appears in the basis.

   double volumetol

        A column is pivoted into the basis when it increases the absolute value
        of the determinant of the basis matrix by more than a factor volumetol.
        This parameter must be >= 1.0. In pratcice typical values are 2.0, 10.0
        or even 100.0. The closer the tolerances to 1.0, the more basis changes
        will usually be necessary to find a maximum volume basis for this
        tolerance (using repeated calls to basiclu_obj_maxvolume(), see below).

    lu_int *p_nupdate

        On return *p_nupdate holds the number of basis updates performed. When
        this is zero and BASICLU_OK is returned, then the volume of the initial
        basis is locally (within one basis change) maximum up to a factor
        volumetol. To find such a basis, basiclu_obj_maxvolume() must be called
        repeatedly starting from an arbitrary basis until *p_nupdate is zero.
        This will happen eventually because each basis update strictly increases
        the volume of the basis matrix. Hence a basis cannot repeat.

        p_nupdate can be NULL, in which case it is not accessed. This is not an
        error condition. The number of updates performed can be obtained as the
        increment to obj->xstore[BASICLU_NUPDATE_TOTAL] caused by the call to
        basiclu_obj_maxvolume().
*/

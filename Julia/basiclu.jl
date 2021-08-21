#
# Julia interface to BASICLU
#

module basiclu

using LinearAlgebra
using LinearAlgebra: checksquare
using Printf
using SparseArrays
using SparseArrays: getcolptr
using Libdl

if haskey(ENV, "JULIA_BASICLU_LIBRARY_PATH")
  const libbasiclu = joinpath(ENV["JULIA_BASICLU_LIBRARY_PATH"], "libbasiclu.$dlext")
else
  using basiclu_jll
end

# comment out if the MAT module is installed
# include("test.jl")

# status codes
const BASICLU_OK = 0
const BASICLU_WARNING_singular_matrix = 2
const BASICLU_ERROR_invalid_store = -1
const BASICLU_ERROR_invalid_call = -2
const BASICLU_ERROR_argument_missing = -3
const BASICLU_ERROR_invalid_argument = -4
const BASICLU_ERROR_maximum_updates = -5
const BASICLU_ERROR_singular_update = -6
const BASICLU_ERROR_invalid_object = -8
const BASICLU_ERROR_out_of_memory = -9

# user parameters in xstore (zero based)
const BASICLU_DROP_TOLERANCE = 4
const BASICLU_ABS_PIVOT_TOLERANCE = 5
const BASICLU_REL_PIVOT_TOLERANCE = 6
const BASICLU_BIAS_NONZEROS = 7
const BASICLU_MAXN_SEARCH_PIVOT = 8
const BASICLU_PAD = 9
const BASICLU_STRETCH = 10
const BASICLU_COMPRESSION_THRESHOLD = 11
const BASICLU_SPARSE_THRESHOLD = 12
const BASICLU_REMOVE_COLUMNS = 13
const BASICLU_SEARCH_ROWS = 14

# user readable from xstore (zero based)
const BASICLU_DIM = 64
const BASICLU_STATUS = 65
const BASICLU_NUPDATE = 70
const BASICLU_NFORREST = 71
const BASICLU_NFACTORIZE = 72
const BASICLU_NUPDATE_TOTAL = 73
const BASICLU_NFORREST_TOTAL = 74
const BASICLU_NSYMPERM_TOTAL = 75
const BASICLU_LNZ = 76
const BASICLU_UNZ = 77
const BASICLU_RNZ = 78
const BASICLU_MIN_PIVOT = 79
const BASICLU_MAX_PIVOT = 80
const BASICLU_UPDATE_COST = 81
const BASICLU_TIME_FACTORIZE = 82
const BASICLU_TIME_SOLVE = 83
const BASICLU_TIME_UPDATE = 84
const BASICLU_TIME_FACTORIZE_TOTAL = 85
const BASICLU_TIME_SOLVE_TOTAL = 86
const BASICLU_TIME_UPDATE_TOTAL = 87
const BASICLU_LFLOPS = 88
const BASICLU_UFLOPS = 89
const BASICLU_RFLOPS = 90
const BASICLU_CONDEST_L = 91
const BASICLU_CONDEST_U = 92
const BASICLU_MAX_ETA = 93
const BASICLU_NORM_L = 94
const BASICLU_NORM_U = 95
const BASICLU_NORMEST_LINV = 96
const BASICLU_NORMEST_UINV = 97
const BASICLU_MATRIX_ONENORM = 98
const BASICLU_MATRIX_INFNORM = 99
const BASICLU_RESIDUAL_TEST = 111
const BASICLU_MATRIX_NZ = 100
const BASICLU_RANK = 101
const BASICLU_BUMP_SIZE = 102
const BASICLU_BUMP_NZ = 103
const BASICLU_NSEARCH_PIVOT = 104
const BASICLU_NEXPAND = 105
const BASICLU_NGARBAGE = 106
const BASICLU_FACTOR_FLOPS = 107
const BASICLU_TIME_SINGLETONS = 108
const BASICLU_TIME_SEARCH_PIVOT = 109
const BASICLU_TIME_ELIM_PIVOT = 110
const BASICLU_PIVOT_ERROR = 120

mutable struct Param
    # Factorization parameters
    absPivotTol::Float64
    relPivotTol::Float64
    dropTol::Float64
    biasNz::Int64
    maxSearch::Int64
    searchRows::Int64
    removeCols::Int64
    memPad::Int64
    memStretch::Float64
    memCompress::Float64
    # Solve parameters
    sparseThres::Float64
    function Param()
        param = new()
        for f in fieldnames(Param)
            setfield!(param, f, convert(fieldtype(Param, f), 0))
        end
        param
    end
end

mutable struct Info
    dim::Int64
    status::Int64
    # Info from last factorization
    nPivot::Int64
    nMatElem::Int64
    nBumpCol::Int64
    nBumpElem::Int64
    nSearchPivot::Int64
    nExpand::Int64
    nGarbage::Int64
    nFactorFlop::Int64
    timeFactorize::Float64
    timeSingletons::Float64
    timeSearchPivot::Float64
    timeElimPivot::Float64
    residualTest::Float64
    onenormMatrix::Float64
    infnormMatrix::Float64
    normL::Float64
    normU::Float64
    normestLinv::Float64
    normestUinv::Float64
    condestL::Float64
    condestU::Float64
    # Info from solves and updates since last factorization
    nUpdate::Int64
    nForrest::Int64
    nElemL::Int64
    nElemU::Int64
    nElemR::Int64
    nFlopL::Int64
    nFlopU::Int64
    nFlopR::Int64
    updateCost::Float64
    minPivot::Float64
    maxPivot::Float64
    maxEta::Float64
    pivotError::Float64
    timeSolve::Float64
    timeUpdate::Float64
    # Accumulated info since object was created
    nFactorize::Int64
    nUpdateTotal::Int64
    nForrestTotal::Int64
    nSympermTotal::Int64
    timeFactorizeTotal::Float64
    timeSolveTotal::Float64
    timeUpdateTotal::Float64
    function Info()
        info = new()
        for f in fieldnames(Info)
            setfield!(info, f, convert(fieldtype(Info, f), 0))
        end
        info
    end
end

function paramnumber(name::Symbol)
    if name == :absPivotTol
        return BASICLU_ABS_PIVOT_TOLERANCE
    elseif name == :relPivotTol
        return BASICLU_REL_PIVOT_TOLERANCE
    elseif name == :dropTol
        return BASICLU_DROP_TOLERANCE
    elseif name == :biasNz
        return BASICLU_BIAS_NONZEROS
    elseif name == :maxSearch
        return BASICLU_MAXN_SEARCH_PIVOT
    elseif name == :searchRows
        return BASICLU_SEARCH_ROWS
    elseif name == :removeCols
        return BASICLU_REMOVE_COLUMNS
    elseif name == :memPad
        return BASICLU_PAD
    elseif name == :memStretch
        return BASICLU_STRETCH
    elseif name == :memCompress
        return BASICLU_COMPRESSION_THRESHOLD
    elseif name == :sparseThres
        return BASICLU_SPARSE_THRESHOLD
    end
    throw(ErrorException("unknown parameter name"))
end

function infonumber(name::Symbol)
    if name == :dim
        return BASICLU_DIM
    elseif name == :status
        return BASICLU_STATUS
    elseif name == :nPivot
        return BASICLU_RANK
    elseif name == :nMatElem
        return BASICLU_MATRIX_NZ
    elseif name == :nBumpCol
        return BASICLU_BUMP_SIZE
    elseif name == :nBumpElem
        return BASICLU_BUMP_NZ
    elseif name == :nSearchPivot
        return BASICLU_NSEARCH_PIVOT
    elseif name == :nExpand
        return BASICLU_NEXPAND
    elseif name == :nGarbage
        return BASICLU_NGARBAGE
    elseif name == :nFactorFlop
        return BASICLU_FACTOR_FLOPS
    elseif name == :timeFactorize
        return BASICLU_TIME_FACTORIZE
    elseif name == :timeSingletons
        return BASICLU_TIME_SINGLETONS
    elseif name == :timeSearchPivot
        return BASICLU_TIME_SEARCH_PIVOT
    elseif name == :timeElimPivot
        return BASICLU_TIME_ELIM_PIVOT
    elseif name == :residualTest
        return BASICLU_RESIDUAL_TEST
    elseif name == :onenormMatrix
        return BASICLU_MATRIX_ONENORM
    elseif name == :infnormMatrix
        return BASICLU_MATRIX_INFNORM
    elseif name == :normL
        return BASICLU_NORM_L
    elseif name == :normU
        return BASICLU_NORM_U
    elseif name == :normestLinv
        return BASICLU_NORMEST_LINV
    elseif name == :normestUinv
        return BASICLU_NORMEST_UINV
    elseif name == :condestL
        return BASICLU_CONDEST_L
    elseif name == :condestU
        return BASICLU_CONDEST_U
    elseif name == :nUpdate
        return BASICLU_NUPDATE
    elseif name == :nForrest
        return BASICLU_NFORREST
    elseif name == :nElemL
        return BASICLU_LNZ
    elseif name == :nElemU
        return BASICLU_UNZ
    elseif name == :nElemR
        return BASICLU_RNZ
    elseif name == :nFlopL
        return BASICLU_LFLOPS
    elseif name == :nFlopU
        return BASICLU_UFLOPS
    elseif name == :nFlopR
        return BASICLU_RFLOPS
    elseif name == :updateCost
        return BASICLU_UPDATE_COST
    elseif name == :minPivot
        return BASICLU_MIN_PIVOT
    elseif name == :maxPivot
        return BASICLU_MAX_PIVOT
    elseif name == :maxEta
        return BASICLU_MAX_ETA
    elseif name == :pivotError
        return BASICLU_PIVOT_ERROR
    elseif name == :timeSolve
        return BASICLU_TIME_SOLVE
    elseif name == :timeUpdate
        return BASICLU_TIME_UPDATE
    elseif name == :nFactorize
        return BASICLU_NFACTORIZE
    elseif name == :nUpdateTotal
        return BASICLU_NUPDATE_TOTAL
    elseif name == :nForrestTotal
        return BASICLU_NFORREST_TOTAL
    elseif name == :nSympermTotal
        return BASICLU_NSYMPERM_TOTAL
    elseif name == :timeFactorizeTotal
        return BASICLU_TIME_FACTORIZE_TOTAL
    elseif name == :timeSolveTotal
        return BASICLU_TIME_SOLVE_TOTAL
    elseif name == :timeUpdateTotal
        return BASICLU_TIME_UPDATE_TOTAL
    end
    throw(ErrorException("unknown info name"))
end

mutable struct basiclu_object
    istore::Ptr{Int64}
    xstore::Ptr{Float64}
    Li::Ptr{Int64}
    Ui::Ptr{Int64}
    Wi::Ptr{Int64}
    Lx::Ptr{Float64}
    Ux::Ptr{Float64}
    Wx::Ptr{Float64}
    lhs::Ptr{Float64}
    ilhs::Ptr{Int64}
    nzlhs::Int64
    realloc_factor::Float64
    function basiclu_object(dim::Int64)
        obj = new()
        retcode = ccall((:basiclu_obj_initialize, libbasiclu), Int64,
                        (Ptr{basiclu_object}, Int64),
                        pointer_from_objref(obj), dim)
        checkretcode("basiclu_obj_initialize", retcode)
        finalizer(obj) do x
            ccall((:basiclu_obj_free, libbasiclu), Cvoid,
                  (Ptr{basiclu_object},), pointer_from_objref(x))
        end
    end
end

function getdim(obj::basiclu_object)
    ccall((:basiclu_obj_get_dim, libbasiclu), Int64,
          (Ptr{basiclu_object},), pointer_from_objref(obj))
end

function getparam(obj::basiclu_object, name::Symbol)
    val = 0.0
    if obj.xstore != C_NULL
        val = unsafe_load(obj.xstore, paramnumber(name) + 1)
    end
    convert(fieldtype(Param, name), val)
end

function getparam(obj::basiclu_object)
    param = Param()
    for f in fieldnames(Param)
        setfield!(param, f, getparam(obj, f))
    end
    param
end

function setparam!(obj::basiclu_object, name::Symbol, val)
    if obj.xstore != C_NULL
        unsafe_store!(obj.xstore, val, paramnumber(name) + 1)
    end
end

function setparam!(obj::basiclu_object, param::Param)
    for f in fieldnames(Param)
        setparam!(obj, f, getfield(param, f))
    end
end

function getinfo(obj::basiclu_object, name::Symbol)
    val = 0.0
    if obj.xstore != C_NULL
        val = unsafe_load(obj.xstore, infonumber(name) + 1)
    end
    convert(fieldtype(Info, name), val)
end

function getinfo(obj::basiclu_object)
    info = Info()
    for f in fieldnames(Info)
        setfield!(info, f, getinfo(obj, f))
    end
    info
end

"""
    factorize(obj::basiclu_object, B::SparseMatrixCSC{Float64, Int64}; check::Bool=true)

Factorize sparse matrix `B`, which must be square and have the dimension for
which `obj` was created.

The factorization method stops when all elements of the active submatrix are
below the absolute pivot tolerance. In this case the `L` and `U` matrix are
padded with columns of the identity matrix, yielding the factorization of a
matrix `B*` which equals `B` except that the dependent columns are replaced by
columns of the identity matrix. The number of actual pivot steps performed can
be obtained with `getinfo(obj, :nPivot)`.

When `check = true`, an error is thrown if the number of pivot steps is less
than the dimension of `B`.
"""
function factorize(obj::basiclu_object, B::SparseMatrixCSC{Float64, Int64}; check::Bool=true)
    dim = getdim(obj)
    m = checksquare(B)
    if dim != m
        throw(DimensionMismatch("matrix dimension does not match basiclu object"))
    end
    if dim != 0
        Bp = getcolptr(B) .- 1
        Bi = rowvals(B) .- 1
        Bx = nonzeros(B)        # don't need a copy
        retcode = ccall((:basiclu_obj_factorize, libbasiclu), Int64,
                        (Ptr{basiclu_object}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
                        pointer_from_objref(obj), Bp, pointer(Bp, 2), Bi, Bx)
        if retcode == BASICLU_WARNING_singular_matrix && !check
            retcode = BASICLU_OK
        end
        checkretcode("basiclu_obj_factorize", retcode)
    end
end

"""
    getfactors(obj::basiclu_object) -> L::SparseMatrixCSC{Float64, Int64},
                                       U::SparseMatrixCSC{Float64, Int64},
                                       p::Vector{Int64},
                                       q::Vector{Int64}

Extract LU factors after fresh factorization. `L` is unit lower triangular, `U`
is upper triangular and `p` and `q` are permutation vectors such that (ignoring
round-off errors) `B[p,q]=LU` when matrix `B` was factorizied.
"""
function getfactors(obj::basiclu_object)
    if getinfo(obj, :nUpdate) != 0
        throw(ErrorException("cannot extract LU factors after factorization has been updated"))
    end
    dim = getdim(obj)
    if dim == 0
        L = spzeros(0, 0)
        U = spzeros(0, 0)
        rowperm = zeros(Int64, 0)
        colperm = zeros(Int64, 0)
        return L, U, rowperm, colperm
    end
    info = getinfo(obj)
    rowperm = Vector{Int64}(undef, dim)
    colperm = Vector{Int64}(undef, dim)
    Lcolptr = Vector{Int64}(undef, dim + 1)
    Lrowidx = Vector{Int64}(undef, info.nElemL + dim)
    Lvalue = Vector{Float64}(undef, info.nElemL + dim)
    Ucolptr = Vector{Int64}(undef, dim + 1)
    Urowidx = Vector{Int64}(undef, info.nElemU + dim)
    Uvalue = Vector{Float64}(undef, info.nElemU + dim)
    retcode = ccall((:basiclu_obj_get_factors, libbasiclu), Int64,
                    (Ptr{basiclu_object}, Ptr{Int64}, Ptr{Int64},
                     Ptr{Int64}, Ptr{Int64}, Ptr{Float64},
                     Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
                    pointer_from_objref(obj), rowperm, colperm,
                    Lcolptr, Lrowidx, Lvalue,
                    Ucolptr, Urowidx, Uvalue)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no fresh factorization available"))
    end
    checkretcode("basiclu_get_factors", retcode)
    L = SparseMatrixCSC{Float64, Int64}(dim, dim, Lcolptr .+ 1, Lrowidx .+ 1, Lvalue)
    U = SparseMatrixCSC{Float64, Int64}(dim, dim, Ucolptr .+ 1, Urowidx .+ 1, Uvalue)
    return L, U, rowperm .+ 1, colperm .+ 1
end

"""
    solve(obj::basiclu_object, rhs::Vector{Float64}, trans::AbstractChar) -> Vector{Float64}

Solve linear system with dense right-hand side. `rhs` is not modified. `trans`
must be 'T' for transposed solve or 'N' for forward solve.
"""
function solve(obj::basiclu_object, rhs::Vector{Float64}, trans::AbstractChar)
    lhs = copy(rhs)
    solve!(obj, lhs, trans)
end

"""
    solve!(obj::basiclu_object, rhs::Vector{Float64}, trans::AbstractChar) -> Vector{Float64}

Solve linear system with dense right-hand side. Solution overwrites `rhs`.
`trans` must be 'T' for transposed solve or 'N' for forward solve.
"""
function solve!(obj::basiclu_object, rhs::Vector{Float64}, trans::AbstractChar)
    checktrans(trans)
    dim = getdim(obj)
    if length(rhs) != dim
        throw(DimensionMismatch("dimension of right-hand side does not match basiclu object"))
    end
    if dim != 0
        retcode = ccall((:basiclu_obj_solve_dense, libbasiclu), Int64,
                        (Ptr{basiclu_object},  Ptr{Float64}, Ptr{Float64}, Cchar),
                        pointer_from_objref(obj), rhs, rhs, trans)
        if retcode == BASICLU_ERROR_invalid_call
            throw(ErrorException("no factorization available"))
        end
        checkretcode("basiclu_obj_solve_dense", retcode)
    end
    rhs
end

"""
    solve(obj::basiclu_object, rhs::SparseVector{Float64, Int64}, trans::AbstractChar) -> SparseVector{Float64, Int64}

Solve linear system with sparse right-hand side. `rhs` is not modified. `trans`
must be 'T' for transposed solve or 'N' for forward solve.
"""
function solve(obj::basiclu_object, rhs::SparseVector{Float64, Int64}, trans::AbstractChar)
    lhs = copy(rhs)
    solve!(obj, lhs, trans)
end

"""
    solve!(obj::basiclu_object, rhs::SparseVector{Float64, Int64}, trans::AbstractChar) -> SparseVector{Float64, Int64}

Solve linear system with sparse right-hand side. Solution overwrites `rhs`.
`trans` must be 'T' for transposed solve or 'N' for forward solve.
"""
function solve!(obj::basiclu_object, rhs::SparseVector{Float64, Int64}, trans::AbstractChar)
    checktrans(trans)
    dim = getdim(obj)
    if length(rhs) != dim
        throw(DimensionMismatch("dimension of right-hand side does not match basiclu object"))
    end
    if dim == 0
        return spzeros(dim)
    end
    retcode = ccall((:basiclu_obj_solve_sparse, libbasiclu), Int64,
                    (Ptr{basiclu_object},  Int64, Ptr{Int64}, Ptr{Float64}, Cchar),
                    pointer_from_objref(obj), nnz(rhs), rowvals(rhs) .- 1, nonzeros(rhs), trans)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no factorization available"))
    end
    checkretcode("basiclu_obj_solve_sparse", retcode)
    gathersol!(obj, rhs)
end

"""
    solve_for_update(obj::basiclu_object, rhs::SparseVector{Float64, Int64}; getsol::Bool=false) -> SparseVector{Float64, Int64}

Solve forward system in preparation to update the factorization. `rhs` holds the
column to be inserted into the factorized matrix in the next call to
[`update`](@ref). When `getsol = true`, then the solution from the forward
solve with right-hand side `rhs` is returned. Otherwise only the update is
prepared.
"""
function solve_for_update(obj::basiclu_object, rhs::SparseVector{Float64, Int64}; getsol::Bool=false)
    dim = getdim(obj)
    if length(rhs) != dim
        throw(DimensionMismatch("dimension of right-hand side does not match basiclu object"))
    end
    if dim == 0
        return getsol ? spzeros(dim) : nothing
    end
    retcode = ccall((:basiclu_obj_solve_for_update, libbasiclu), Int64,
                    (Ptr{basiclu_object}, Int64, Ptr{Int64}, Ptr{Float64}, Cchar, Int64),
                    pointer_from_objref(obj), nnz(rhs), rowvals(rhs) .- 1, nonzeros(rhs), 'N', getsol ? 1 : 0)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no factorization available"))
    end
    if retcode == BASICLU_ERROR_maximum_updates
        throw(ErrorException("maximum number of updates reached"))
    end
    checkretcode("basiclu_obj_solve_for_update", retcode)
    return getsol ? gathersol(obj) : nothing
end

"""
    solve_for_update(obj::basiclu_object, pos::Int64; getsol::Bool=false) -> SparseVector{Float64, Int64}

Solve transposed system in preparation to update the factorization. `pos` holds
the column index of the factorized matrix to be replaced in the next call to
[`update`](@ref). When `getsol = true`, then the solution from the transposed
solve with a unit vector as right-hand side is returned. Otherwise only the
update is prepared.
"""
function solve_for_update(obj::basiclu_object, pos::Int64; getsol::Bool=false)
    dim = getdim(obj)
    if pos < 1 || pos > dim
        throw(DimensionMismatch("column index outside dimension of basiclu object"))
    end
    if dim == 0
        return getsol ? spzeros(dim) : nothing
    end
    retcode = ccall((:basiclu_obj_solve_for_update, libbasiclu), Int64,
                    (Ptr{basiclu_object}, Int64, Ptr{Int64}, Ptr{Float64}, Cchar, Int64),
                    pointer_from_objref(obj), 0, Ref{Int64}(pos-1), C_NULL, 'T', getsol ? 1 : 0)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no factorization available"))
    end
    if retcode == BASICLU_ERROR_maximum_updates
        throw(ErrorException("maximum number of updates reached"))
    end
    checkretcode("basiclu_obj_solve_for_update", retcode)
    return getsol ? gathersol(obj) : nothing
end

"""
    update(obj::basiclu_object, pivot::Float64) -> Float64

Update the factorization after a column modification. The column position and
the new column must have been set in previous calls to
[`solve_for_update`](@ref).

`pivot` is the pivot element corresponding to the update operation; i.e. when
column `j` of `B` is to be replaced by vector `v`, then `pivot = (B\\v)[j]`. The
absolute difference between `pivot` and a recomputed version can be obtained
with `getinfo(obj, :pivotError)`; this is also the return value. A pivot error
larger than 1e-8, say, indicates numerical instability and suggests
refactorization.

An error is thrown when the recomputed pivot element is below the absolute pivot
tolerance. In this case no update is performed and the old factorization remains
valid.
"""
function update(obj::basiclu_object, pivot::Float64)
    dim = getdim(obj)
    if dim != 0
        retcode = ccall((:basiclu_obj_update, libbasiclu), Int64,
                        (Ptr{basiclu_object},  Float64),
                        pointer_from_objref(obj), pivot)
        if retcode == BASICLU_ERROR_invalid_call
            throw(ErrorException("not prepared for update"))
        end
        checkretcode("basiclu_update", retcode)
    end
    getinfo(obj, :pivotError)
end

"""
    maxvolume(obj::basiclu_object, A::SparseMatrixCSC{Float64, Int64}, basis::Vector{Int64}, volumetol::Float64=2.0) -> Int64

Given an initial basis such that `A[:,basis]` is square and nonsingular, make
one pass over the nonbasic columns of `A` and pivot each column into the basis
when it increases the absolute value of the determinant of the basis matrix by
more than a factor `volumetol`. On return `basis` has been updated. Return the
number of basis updates performed.
"""
function maxvolume(obj::basiclu_object, A::SparseMatrixCSC{Float64, Int64}, basis::Vector{Int64},
                   volumetol::Float64=2.0)
    m, n = size(A)
    if m != getdim(obj)
        throw(DimensionMismatch("number of rows of A does not match dimension of basiclu object"))
    end
    if length(basis) != getdim(obj)
        throw(DimensionMismatch("basis has incorrect length"))
    end
    isbasic = zeros(Int64, n)
    isbasic[basis] .= 1
    if sum(isbasic) != m
        throw(ErrorException("duplicate index in basis"))
    end
    cbasis = basis .- 1
    Ap = getcolptr(A) .- 1
    Ai = rowvals(A) .- 1
    Ax = nonzeros(A)            # don't need a copy
    p_nupdate = Ref{Int64}(0)
    retcode = ccall((:basiclu_obj_maxvolume, libbasiclu), Int64,
                    (Ptr{basiclu_object}, Int64, Ptr{Int64}, Ptr{Int64}, Ptr{Float64},
                     Ptr{Int64}, Ptr{Int64}, Float64, Ptr{Int64}),
                    pointer_from_objref(obj), n, Ap, Ai, Ax, cbasis, isbasic, volumetol, p_nupdate)
    basis[:] = cbasis .+ 1
    checkretcode("basiclu_obj_maxvolume", retcode)
    return p_nupdate[]
end

"""
    maxvolbasis(A::SparseMatrixCSC{Float64, Int64}; lindeptol::Float64=1e-8,
                volumetol::Float64=2.0, maxpass::Int64=2, verbose::Bool=true) -> Vector{Int64}, basiclu_object 

Find a set of column indices for the matrix `AI = [A I]` such that `AI[:,basis]`
is square and nonsingular and the number of slack columns in the basis is
minimum (this is the row rank deficiency of `A`). Return the vector of column
indices of `AI` which form the basis matrix and a `basiclu_object` which holds a
factorization of the basis matrix.

Method: Scale the slack columns of `AI` by `lindeptol` and try to find a maximum
volume basis for this matrix by making at most `maxpass` calls to
[`maxvolume`](@ref). If `verbose` is true, then print the number of basis
updates after each call.
"""
function maxvolbasis(A::SparseMatrixCSC{Float64, Int64}; lindeptol::Float64=1e-8,
                     volumetol::Float64=2.0, maxpass::Int64=2,
                     verbose::Bool=true)
    m, n = size(A)
    rowmax = maximum(abs.(A), dims=2)[:]
    rowmax = max.(rowmax, 1.)
    AI = [A spdiagm(0 => lindeptol.*rowmax)]
    basis = collect(n+1:n+m)
    obj = basiclu_object(m)
    for pass = 1:maxpass
        nupdate = maxvolume(obj, AI, basis, volumetol)
        if verbose
            @printf("pass %d: %d updates\n", pass, nupdate)
        end
        if nupdate == 0 break; end
    end
    return basis, obj
end

function checktrans(trans::AbstractChar)
    if !(trans == 'N' || trans == 'T')
        throw(ArgumentError("trans argument must be 'N' (no transpose) or 'T' (transpose), got '$trans'"))
    end
    trans
end

function checkretcode(funcname::String, retcode::Int64)
    if retcode == BASICLU_ERROR_out_of_memory
        throw(OutOfMemoryError())
    elseif retcode == BASICLU_WARNING_singular_matrix
        msg = @sprintf("%s() encountered singular matrix", funcname)
        throw(ErrorException(msg))
    elseif retcode == BASICLU_ERROR_singular_update
        msg = @sprintf("%s() encountered singular update", funcname)
        throw(ErrorException(msg))
    elseif retcode != BASICLU_OK
        msg = @sprintf("%s() failed with return code %d", funcname, retcode)
        throw(ErrorException(msg))
    end
end

function gathersol(obj::basiclu_object)
    dim = getdim(obj)
    nzind = Vector{Int64}(undef, obj.nzlhs)
    nzval = Vector{Float64}(undef, obj.nzlhs)
    lhs = SparseVector{Float64, Int64}(dim, nzind, nzval)
    gathersol!(obj, lhs)
end

function gathersol!(obj::basiclu_object, lhs::SparseVector{Float64, Int64})
    dim = getdim(obj)
    @assert length(lhs) == dim
    if nnz(lhs) != obj.nzlhs
        resize!(lhs.nzind, obj.nzlhs)
        resize!(lhs.nzval, obj.nzlhs)
    end
    for i = 1:obj.nzlhs
        idx = unsafe_load(obj.ilhs, i) + 1
        @assert idx >= 1 && idx <= dim
        lhs.nzind[i] = idx
        lhs.nzval[i] = unsafe_load(obj.lhs, idx)
    end
    p = sortperm(lhs.nzind)
    lhs.nzind[:] = lhs.nzind[p]
    lhs.nzval[:] = lhs.nzval[p]
    lhs
end

end

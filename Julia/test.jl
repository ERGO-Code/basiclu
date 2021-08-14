using MAT

function readmat(file, names...)
    n = length(names)
    objects = Array{Any}(undef,n)
    f1 = matopen(file)
    for k = 1:n
        objects[k] = read(f1, names[k])
    end
    close(f1)
    return tuple(objects...)
end

"""
    test_small()

Test special cases.
"""
function test_small()
    # dimension 0
    B = spzeros(0, 0)
    blu = basiclu_object(0)
    factorize(blu, B)
    rhs = zeros(0)
    lhs = solve(blu, rhs, 'N')
    lhs = solve(blu, rhs, 'T')
    rhs = spzeros(0)
    lhs = solve(blu, rhs, 'N')
    lhs = solve(blu, rhs, 'T')
    lhs = solve_for_update(blu, rhs, getsol=true)

    # dimension 1
    B = sparse(1.0I, 1, 1)
    blu = basiclu_object(1)
    factorize(blu, B)
    rhs = [2.0]
    lhs = solve(blu, rhs, 'N')
    @assert norm(lhs-rhs, Inf) == 0.0
    lhs = solve(blu, rhs, 'T')
    @assert norm(lhs-rhs, Inf) == 0.0
    rhs = SparseVector{Float64,Int64}(1, [1], [2.0])
    lhs = solve(blu, rhs, 'N')
    @assert norm(lhs-rhs, Inf) == 0.0
    lhs = solve(blu, rhs, 'T')
    @assert norm(lhs-rhs, Inf) == 0.0
    lhs = solve_for_update(blu, rhs, getsol=true)
    @assert norm(lhs-rhs, Inf) == 0.0

    # dimension 1 singular
    B = spzeros(1, 1)
    blu = basiclu_object(1)
    factorize(blu, B, check=false)
    info = getinfo(blu)
    @assert info.nPivot == 0
    # basiclu has inserted a slack column for the zero pivot
    rhs = [2.0]
    lhs = solve(blu, rhs, 'N')
    @assert norm(lhs-rhs, Inf) == 0.0
    lhs = solve(blu, rhs, 'T')
    @assert norm(lhs-rhs, Inf) == 0.0
    rhs = SparseVector{Float64,Int64}(1, [1], [2.0])
    lhs = solve(blu, rhs, 'N')
    @assert norm(lhs-rhs, Inf) == 0.0
    lhs = solve(blu, rhs, 'T')
    @assert norm(lhs-rhs, Inf) == 0.0
    lhs = solve_for_update(blu, rhs, getsol=true)
    @assert norm(lhs-rhs, Inf) == 0.0
end

"""
    test_factorize(testdir::String, trans::Bool=false)

For all `*.mat` files in `testdir` read matrix `B` and factorize it. Monitor
residual of factorization and forward/backward solves.

`trans` specifies whether `B` or its transposed is factorized.
"""
function test_factorize(testdir::String, trans::Bool=false)
    files = readdir(testdir)
    for f in files
        if length(f) < 4 || f[end-3:end] != ".mat"
            continue
        end
        @printf(" %-24s", f)
        B, = readmat(string(testdir, f), "B")
        Bt = sparse(transpose(B))
        if trans
            (B,Bt) = (Bt,B)
        end
        m = size(B,1)
        blu = basiclu_object(m)
        factorize(blu, B)
        L,U,p,q = getfactors(blu)
        res = norm(L*U-B[p,q], Inf)

        # test dense solve
        rhs = ones(m)
        lhs = solve(blu, rhs, 'N')
        res = max(res, norm(B*lhs-rhs, Inf))
        lhs = solve(blu, rhs, 'T')
        res = max(res, norm(Bt*lhs-rhs, Inf))

        # test sparse solve
        rhs = SparseVector{Float64,Int64}(m, [1;m], [1.0;1.0])
        lhs = solve(blu, rhs, 'N')
        res = max(res, norm(B*lhs-rhs, Inf))
        lhs = solve(blu, rhs, 'T')
        res = max(res, norm(Bt*lhs-rhs, Inf))

        # test solve for update
        lhs = solve_for_update(blu, rhs, getsol=true)
        res = max(res, norm(B*lhs-rhs, Inf))
        lhs = solve_for_update(blu, m-1, getsol=true)
        rhs = zeros(m); rhs[m-1] = 1
        res = max(res, norm(Bt*lhs-rhs, Inf))
        @printf("%.2e\n", res)
    end
    nothing
end

"""
    test_update(testdir::String)

For all `*.mat` files in `testdir` read matrix `A` and the pivot sequence.
Starting from the slack basis apply pivot operations. Monitor residuals to
forward/backward solves after each 100 updates.
"""
function test_update(testdir::String)
    files = readdir(testdir)
    for f in files
        if length(f) < 4 || f[end-3:end] != ".mat"
            continue
        end
        @printf(" %-24s", f)
        A,invar,outvar = readmat(string(testdir, f), "A", "invar", "outvar")
        invar = invar[:]
        outvar = outvar[:]
        m,n = size(A)
        A1 = [A sparse(1.0I, m, m)]
        basis = Vector{Int64}(undef,m)
        basis[:] = collect(1:m) .+ n # slack basis
        res,nfactor,nforrest,nperm = test_update(A1, basis, invar, outvar)
        @printf("%.2e      %5d factor, %6d forrest, %6d perm\n",
                res, nfactor, nforrest, nperm)
    end
    nothing
end

function test_update(A::SparseMatrixCSC{Float64,Int64}, basis::Vector{Int64}, invar::Vector{Int64},
                     outvar::Vector{Int64})
    m,n = size(A)
    niter = length(invar)
    map2basis = zeros(Int64, n)
    map2basis[basis] = 1:m
    blu = basiclu_object(m)
    param = getparam(blu)
    #param.relPivotTol = 0.5
    setparam!(blu, param)
    factorize(blu, A[:,basis])
    res = 0.0
    for k = 1:niter
        if invar[k] == outvar[k] # bound flip
            continue
        end
        j = invar[k]
        p = map2basis[outvar[k]]
        @assert p > 0
        lhs = solve_for_update(blu, A[:,j], getsol=true)
        piv = lhs[p]
        solve_for_update(blu, p)
        piverr = update(blu, piv)
        basis[p] = j
        map2basis[j] = p
        map2basis[outvar[k]] = 0
        cost = getinfo(blu, :updateCost)
        if piverr > 1e-10 || cost > 1.0
            factorize(blu, A[:,basis])
        end
        if k%100 == 0
            B = A[:,basis]
            rhs = SparseVector{Float64,Int64}(m, [1;m], [1.0;1.0])
            lhs = solve(blu, rhs, 'N')
            res = max(res, norm(B*lhs-rhs, Inf))
            lhs = solve(blu, rhs, 'T')
            res = max(res, norm(B'*lhs-rhs, Inf))
        end
    end
    info = getinfo(blu)
    return res, info.nFactorize, info.nForrestTotal, info.nUpdateTotal - info.nForrestTotal
end

"""
    test_maxvolume(testdir::String)

For all `*.mat` files in `testdir` read matrix `A` and call
[`maxvolbasis`](@ref) to construct a basis.
"""
function test_maxvolume(testdir::String)
    files = readdir(testdir)
    for f in files
        if length(f) < 4 || f[end-3:end] != ".mat"
            continue
        end
        @printf("%s\n", f)
        A, = readmat(string(testdir, f), "A")
        m, n = size(A)
        basis, obj = basiclu.maxvolbasis(A, maxpass=5)
        AI = [A I]
        B = AI[:, basis]
        # refactorize to check that basis is nonsingular
        setparam!(obj, :relPivotTol, 0.3)
        factorize(obj, B)
        info = basiclu.getinfo(obj)
        nslack = sum(basis.>n)
        @printf("%d slack, condestL = %.2e, condestU = %.2e\n\n",
                nslack, info.condestL, info.condestU)
    end
    nothing
end

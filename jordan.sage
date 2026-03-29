from sage.all import *

def jordanhermi_unr(A, alg_data):
    """
    Computes the Jordan decomposition of a Hermitian matrix over a general 
    local algebraic number field for the unramified (or inert/ramified non-split) case.
    
    INPUT:
        A        : matrix; The Hermitian matrix to be decomposed.
        alg_data : dict; Dictionary containing local algebraic structures.
                        
    OUTPUT:
        tuple: (newbasis, gram_mat, dims)
    """
    dim = A.nrows()
    d = 0
    dims = []
    B = copy(A)
    newbasis = Matrix(alg_data['E'], dim, dim, 1)

    vals_B = [alg_data['val_func'](x) for x in B.list() if x != 0]
    i = min(vals_B) if vals_B else 0

    while d < dim:
        while Matrix(alg_data['residue_field'], B / (alg_data['unif']**i)) == 0:
            i += 1
            dims.append(0)

        Bmod = Matrix(alg_data['residue_field'], B / (alg_data['unif']**i))
        rank = Bmod.rank()
        dims.append(rank)

        kernbase = Bmod.right_kernel().echelonized_basis()
        indices = []
        k = 0
        for v in kernbase:
            l = k + 1
            while v[l-1] == 0:
                l += 1
            indices.extend(range(k+1, l))
            k = l
        indices.extend(range(k+1, dim - d + 1))
        P = Matrix(alg_data['OE'], dim - d, dim - d)
        k = 1
        for j in indices:
            P[j-1, k-1] = 1
            k += 1

        kernbase_oe = (B * P).integer_kernel(alg_data['OE']).basis()
        for v in kernbase_oe:
            l = 0
            for x in v:
                P[l, k-1] = alg_data['conj'](x)
                l += 1
            k += 1

        if d + rank < dim:
            Ptransconj = Matrix(alg_data['OE'], [[alg_data['conj'](P[b, a]) for b in range(dim - d)] for a in range(dim - d)])
            B = (Ptransconj * B * P).submatrix(rank, rank)

        newbasis = newbasis * block_diagonal_matrix(diagonal_matrix([1]*d), P, subdivide=False)
        d += rank
        i += 1

    newbasistransconj = Matrix(alg_data['E'], [[alg_data['conj'](newbasis[b, a]) for b in range(dim)] for a in range(dim)])
    return newbasis, newbasistransconj * A * newbasis, dims


def jordanhermi_split(A, alg_data, t):
    """
    Computes the Jordan decomposition of a Hermitian matrix over a general 
    local algebraic number field for the split case.
    
    INPUT:
        A        : matrix; The Hermitian matrix to be decomposed.
        alg_data : dict; Dictionary containing local algebraic structures.
        t        : element; A local idempotent approximator.
                        
    OUTPUT:
        tuple: (newbasis, gram_mat, dims)
    """
    tconj = alg_data['OE'](alg_data['conj'](t))
    
    dim = A.nrows()
    d = 0
    dims = []
    B = copy(A)
    newbasis = Matrix(alg_data['E'], dim, dim, 1)

    vals_B = [alg_data['val_func'](x) for x in B.list() if x != 0]
    i = min(vals_B) if vals_B else 0

    while d < dim:
        while Matrix(alg_data['residue_field'], B / (alg_data['unif']**i)) == 0:
            i += 1
            dims.append(0)

        Bmod = Matrix(alg_data['residue_field'], B / (alg_data['unif']**i))
        rank = Bmod.rank()
        dims.append(rank)

        rkern = Bmod.right_kernel().echelonized_basis()
        indices = []
        k = 0
        for v in rkern:
            l = k + 1
            while v[l-1] == 0:
                l += 1
            indices.extend(range(k+1, l))
            k = l
        indices.extend(range(k+1, dim - d + 1))
        
        P1 = Matrix(alg_data['OE'], dim - d, dim - d)
        k = 1
        for j in indices:
            P1[j-1, k-1] = 1
            k += 1

        lkern = Bmod.left_kernel().echelonized_basis()
        indices = []
        k = 0
        for v in lkern:
            l = k + 1
            while v[l-1] == 0:
                l += 1
            indices.extend(range(k+1, l))
            k = l
        indices.extend(range(k+1, dim - d + 1))
        
        P2 = Matrix(alg_data['OE'], dim - d, dim - d)
        k = 1
        for j in indices:
            P2[j-1, k-1] = 1
            k += 1

        P = (1-t)*P1 + t*P2
        Ptransconj = (1-tconj)*P1.transpose() + tconj*P2.transpose()

        kernbase_oe = (B * P).integer_kernel(alg_data['OE']).basis()
        for v in kernbase_oe:
            l = 0
            for x in v:
                P[l, k-1] = alg_data['conj'](x)
                Ptransconj[k-1, l] = x
                l += 1
            k += 1

        if d + rank < dim:
            B = (Ptransconj * B * P).submatrix(rank, rank)

        newbasis = newbasis * block_diagonal_matrix(diagonal_matrix([1]*d), P, subdivide=False)
        d += rank
        i += 1

    newbasistransconj = Matrix(alg_data['E'], [[alg_data['conj'](newbasis[b, a]) for b in range(dim)] for a in range(dim)])
    return newbasis, newbasistransconj * A * newbasis, dims


def jordanhermi(A, is_split, alg_data, t=None):
    """
    Unified wrapper function for generalized Jordan decomposition of a Hermitian matrix.
    
    INPUT:
        A        : matrix; The Hermitian matrix to be decomposed.
        is_split : bool; True if the relevant prime ideal splits in the extension E/F.
        alg_data : dict; Dictionary containing local algebraic structures.
        t        : element (Optional); The local idempotent approximator.
        
    OUTPUT:
        tuple: (newbasis, gram_mat, dims)
    """
    if is_split:
        if t is None:
            raise ValueError("Split case requires the idempotent approximator 't'.")
        return jordanhermi_split(A, alg_data, t)
    else:
        return jordanhermi_unr(A, alg_data)



def jordanquadZp(A, p, prec=20):
    """
    Computes the Jordan decomposition of a symmetric matrix (quadratic form) over Z_p.
    
    INPUT:
        A    : matrix; Symmetric matrix over QQ or Z_p.
        p    : integer; Prime number.
        prec : integer; p-adic precision for computation.
        
    OUTPUT:
        newbasis : The basis transformation matrix P.
        gram_mat : The decomposed Gram matrix (P.transpose() * A * P).
        dims     : List of ranks for each p^i component.
    """
    
    dim = A.nrows()
    processed_dim = 0
    p_exponent = 0
    dims = []
    B = copy(A)
    newbasis = identity_matrix(QQ, dim)

    def get_val(x):
        try:
            return x.valuation(p)
        except TypeError:
            return x.valuation()

    while processed_dim < dim:
        while True:
            if all(get_val(val / p**p_exponent) >= 1 for val in B.list() if val != 0):
                p_exponent += 1
                dims.append(0)
            else:
                break

        Bmod = matrix(GF(p), B / p**p_exponent)
        rank = Bmod.rank()
        dims.append(rank)

        kernbase = Bmod.right_kernel().echelonized_basis()
        indices = []
        last_idx = 0
        for v in kernbase:
            l = last_idx + 1
            while v[l-1] == 0:
                l += 1
            indices.extend(range(last_idx + 1, l))
            last_idx = l
        
        indices.extend(range(last_idx + 1, dim - processed_dim + 1))

        P = matrix(ZZ, dim - processed_dim, dim - processed_dim)
        col = 1
        for j in indices:
            P[j-1, col-1] = 1
            col += 1

        try:
            kernbase_zp = (B * P).integer_kernel().basis()
        except AttributeError:
            BP_lifted = matrix(QQ, [[QQ(x.lift()) if hasattr(x, 'lift') else QQ(x) for x in row] for row in (B*P)])
            kernbase_zp = BP_lifted.integer_kernel().basis()

        for v in kernbase_zp:
            for row_idx, val in enumerate(v):
                P[row_idx, col-1] = val
            col += 1

        if processed_dim + rank < dim:
            B = (P.transpose() * B * P).submatrix(rank, rank)

        block_P = block_diagonal_matrix(identity_matrix(QQ, processed_dim), P)
        newbasis = newbasis * block_P
        
        processed_dim += rank
        p_exponent += 1

    gram_mat = newbasis.transpose() * A * newbasis
    return newbasis, gram_mat, dims
from sage.all import *

def jordanhermi_unr(A, E, OE, residue_field, unif, conj, val_func):
    """
    Computes the Jordan decomposition of a Hermitian matrix over a general 
    local algebraic number field for the unramified (or inert/ramified non-split) case.
    
    INPUT:
        A             : matrix; The Hermitian matrix to be decomposed.
        E             : Field; The field to which the components belong (quadratic extension E/F).
        OE            : Ring; The ring of integers of E.
        residue_field : Field; The residue field of the base field F (F_q).
        unif          : element; A uniformizer (\varpi_F) of the base field F.
        conj          : function; The non-trivial involution (conjugation) map of E/F.
        val_func      : function; A function that takes an element of E and returns 
                        its valuation at the corresponding prime ideal.
                        
    OUTPUT:
        tuple: (newbasis, gram_mat, dims)
            - newbasis (matrix): The transformation matrix.
            - gram_mat (matrix): The decomposed Gram matrix (newbasistransconj * A * newbasis).
            - dims (list): A list of dimensions for each Jordan block.
    """
    dim = A.nrows()
    d = 0
    dims = []
    B = copy(A)
    newbasis = Matrix(E, dim, dim, 1)

    # Calculate the minimum valuation using the provided val_func instead of p
    vals_B = [val_func(x) for x in B.list() if x != 0]
    i = min(vals_B) if vals_B else 0

    while d < dim:
        # Scale by unif**i (the uniformizer of the base field) instead of p**i
        while Matrix(residue_field, B / (unif**i)) == 0:
            i += 1
            dims.append(0)

        Bmod = Matrix(residue_field, B / (unif**i))
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

        P = Matrix(OE, dim, dim)
        k = 1
        for j in indices:
            P[j-1, k-1] = 1
            k += 1

        # Note: integer_kernel relies on SageMath's support for rings of integers and algebraic fields
        kernbase_oe = (B * P).integer_kernel(OE).basis()
        for v in kernbase_oe:
            l = 0
            for x in v:
                P[l, k-1] = conj(x)
                l += 1
            k += 1

        if d + rank < dim:
            Ptransconj = Matrix(OE, [[conj(P[b, a]) for b in range(dim)] for a in range(dim)])
            B = (Ptransconj * B * P).submatrix(rank, rank)

        newbasis = newbasis * block_diagonal_matrix(diagonal_matrix([1]*d), P, subdivide=False)
        d += rank
        i += 1

    newbasistransconj = Matrix(E, [[conj(newbasis[b, a]) for b in range(dim)] for a in range(dim)])
    return newbasis, newbasistransconj * A * newbasis, dims


def jordanhermi_split(A, E, OE, residue_field, unif, conj, val_func, t):
    """
    Computes the Jordan decomposition of a Hermitian matrix over a general 
    local algebraic number field for the split case.
    
    INPUT:
        A             : matrix; The Hermitian matrix to be decomposed.
        E             : Field; The field to which the components belong.
        OE            : Ring; The ring of integers of E.
        residue_field : Field; The residue field of the base field F (F_q).
        unif          : element; A uniformizer (\varpi_F) of the base field F.
        conj          : function; The non-trivial involution (conjugation) map of E/F.
        val_func      : function; A function returning the valuation of an element of E.
        t             : element; An element of OE acting as a local idempotent approximator 
                        when the prime ideal of F splits into two prime ideals P1 and P2 in E.
                        (i.e., t = 1 mod P1 and t = 0 mod P2).
                        
    OUTPUT:
        tuple: (newbasis, gram_mat, dims)
            - newbasis (matrix): The transformation matrix.
            - gram_mat (matrix): The decomposed Gram matrix.
            - dims (list): A list of dimensions for each Jordan block.
    """
    tconj = OE(conj(t))
    
    dim = A.nrows()
    d = 0
    dims = []
    B = copy(A)
    newbasis = Matrix(E, dim, dim, 1)

    vals_B = [val_func(x) for x in B.list() if x != 0]
    i = min(vals_B) if vals_B else 0

    while d < dim:
        while Matrix(residue_field, B / (unif**i)) == 0:
            i += 1
            dims.append(0)

        Bmod = Matrix(residue_field, B / (unif**i))
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
        
        P1 = Matrix(OE, dim, dim)
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
        
        P2 = Matrix(OE, dim, dim)
        k = 1
        for j in indices:
            P2[j-1, k-1] = 1
            k += 1

        # Local gluing of bases using the idempotent t
        P = (1-t)*P1 + t*P2
        Ptransconj = (1-tconj)*P1.transpose() + tconj*P2.transpose()

        kernbase_oe = (B * P).integer_kernel(OE).basis()
        for v in kernbase_oe:
            l = 0
            for x in v:
                P[l, k-1] = conj(x)
                Ptransconj[k-1, l] = x
                l += 1
            k += 1

        if d + rank < dim:
            B = (Ptransconj * B * P).submatrix(rank, rank)

        newbasis = newbasis * block_diagonal_matrix(diagonal_matrix([1]*d), P, subdivide=False)
        d += rank
        i += 1

    newbasistransconj = Matrix(E, [[conj(newbasis[b, a]) for b in range(dim)] for a in range(dim)])
    return newbasis, newbasistransconj * A * newbasis, dims


def jordanhermi(A, is_split, E, OE, residue_field, unif, conj, val_func, t=None):
    """
    Unified wrapper function for generalized Jordan decomposition of a Hermitian matrix.
    
    INPUT:
        A             : matrix; The Hermitian matrix to be decomposed.
        is_split      : bool; True if the relevant prime ideal splits in the extension E/F.
        E             : Field; The extension field.
        OE            : Ring; The ring of integers of E.
        residue_field : Field; The residue field of the base field F.
        unif          : element; A uniformizer of the base field F.
        conj          : function; The conjugation map of E/F.
        val_func      : function; The valuation function.
        t             : element (Optional); The local idempotent approximator, required if is_split is True.
        
    OUTPUT:
        tuple: (newbasis, gram_mat, dims) generated by the appropriate decomposition logic.
    """
    if is_split:
        if t is None:
            raise ValueError("Split case requires the idempotent approximator 't'.")
        return jordanhermi_split(A, E, OE, residue_field, unif, conj, val_func, t)
    else:
        return jordanhermi_unr(A, E, OE, residue_field, unif, conj, val_func)
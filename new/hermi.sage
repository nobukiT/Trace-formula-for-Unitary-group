load("hermilatticeclass.sage") 
load("young.sage")
load("utils.sage")

# ========================================================================
# 1. Basic Combinatorics and Linear Algebra
# ========================================================================

def partitions(n, maxpower):
    """
    Generate all partitions of a multiplicity `n` up to a maximum depth `maxpower`.
    
    INPUT:
        n        : integer; the total multiplicity (dimension) to partition
        maxpower : integer; the maximum allowed depth
    OUTPUT:
        A list of lists, where each inner list represents a valid partition type.
    """
    L = [[n]]
    if maxpower == 0:
        return L
    i = n - 1
    while i >= 0:
        for K in partitions(n - i, maxpower - 1):
            L.append([i] + K)
        i -= 1
    return L

def quadQp_from_mat(A, p):
    """
    Generate a local quadratic form (QuadQp) at prime p from a symmetric matrix A.
    
    INPUT:
        A : matrix; a symmetric Gram matrix over the local base field F
        p : integer; the rational prime under consideration

    OUTPUT:
        An instance of the QuadQp class representing the isometry class of the 
        local quadratic space.
    """
    dim = A.nrows()
    disc = A.det()
    if disc == 0:
        raise ValueError(f"Discriminant is zero at prime {p}. Matrix:\n{A}")
    
    Q = QuadraticForm(QQ, 2 * A)
    hasse = Q.hasse_invariant(p)
    return QuadQp(p, dim, disc, 1 if hasse > 0 else -1)

def tracematrix(A, E, F, trace_func, basis):
    """
    Compute the trace matrix of a Hermitian matrix A over the extension E/F.
    
    INPUT:
        A          : matrix; a Hermitian matrix defined over the extension field E
        E          : algebraic field; the quadratic extension field
        F          : algebraic field; the base local field
        trace_func : function; the trace map from E to F
        basis      : list; the basis elements of E over F

    OUTPUT:
        A block matrix over the base field F, representing the trace form of A.
    """
    nrows, ncols = A.dimensions()
    deg = len(basis)
    matlist = []
    
    for a in range(nrows):
        for b in range(ncols):
            entries = [[trace_func(A[a,b] * basis[i] * basis[j]) for j in range(deg)] for i in range(deg)]
            matlist.append(matrix(F, deg, deg, entries))
            
    return block_matrix(nrows, ncols, matlist, subdivide=False)

# ========================================================================
# 2. Combinatorial Types for Hermitian Lattices
# ========================================================================

def enum_hermi_ramoddp_types(d, e, mult, maxpower, mini=0):
    """
    Enumerate the invariants (types) of ramified Hermitian lattices at an odd prime.
    
    INPUT:
        d        : integer; valuation of the different of E/F
        e        : integer; ramification index
        mult     : integer; total multiplicity (dimension)
        maxpower : integer; maximum depth of the Jordan decomposition
        mini     : integer; current minimum depth index (used for recursion, default 0)

    OUTPUT:
        A list of lists, where each inner list contains tuples `[dim, type_mod_2]` 
        representing valid Jordan block invariants.
    """
    if mini > maxpower * e:
        return []
    L = []
    if (mini - d) % 2 != 0:
        i = mult if mult % 2 == 0 else mult + 1
        if mult % 2 == 0:
            L.append([[mult]])
        while i > 1:
            i -= 2
            for J in enum_hermi_ramoddp_types(d, e, mult - i, maxpower, mini + 1):
                L.append([[i]] + J)
    else:
        L = [[[mult, Mod(0, 2)]], [[mult, Mod(1, 2)]]]
        i = mult
        while i > 0:
            i -= 1
            for J in enum_hermi_ramoddp_types(d, e, mult - i, maxpower, mini + 1):
                L.append([[i, Mod(0, 2)]] + J)
                if i > 0:
                    L.append([[i, Mod(1, 2)]] + J)
    return L

def formal_types_hermi_ram2(d, e, mult, maxpower, mini=0):
    """
    Enumerate the invariants (types) of ramified Hermitian lattices at p = 2.
    
    INPUT:
        d        : integer; valuation of the different of E/F
        e        : integer; ramification index
        mult     : integer; total multiplicity (dimension)
        maxpower : integer; maximum depth of the Jordan decomposition
        mini     : integer; current minimum depth index (used for recursion, default 0)

    OUTPUT:
        A list of lists containing tuples `[dim, type_info...]` representing 
        valid Jordan block invariants for the p=2 case.
    """
    if mini > maxpower * e:
        return []
    L = []
    if (mini - d) % 2 != 0:  
        i = mult // 2 if mult % 2 == 0 else (mult + 1) // 2
        if mult % 2 == 0:
            L = [[[mult // 2, Mod(0, 2)]], [[mult // 2, Mod(1, 2)]]]
        while i > 0:
            i -= 1
            for J in formal_types_hermi_ram2(d, e, mult - 2 * i, maxpower, mini + 1):
                L.append([[i, Mod(0, 2)]] + J)
                if i > 0 and (J[0][0] == 0 or J[0][1] == Mod(0, 2)):
                    L.append([[i, Mod(1, 2)]] + J)
    else:  
        L = [[[mult, Mod(1, 2), Mod(0, 2)]], [[mult, Mod(1, 2), Mod(1, 2)]]]
        if mult % 2 == 0:
            L.append([[mult, Mod(0, 2)]])
        i = mult
        while i > 0:
            i -= 1
            for J in formal_types_hermi_ram2(d, e, mult - i, maxpower, mini + 1):
                if i % 2 == 0:
                    L.append([[i, Mod(0, 2)]] + J)
                if i > 0:
                    L.append([[i, Mod(1, 2), Mod(0, 2)]] + J)
                    if J[0][0] == 0 and (len(J) == 1 or J[1][0] == 0 or J[1][1] == 0):
                        L.append([[i, Mod(1, 2), Mod(1, 2)]] + J)
    return L

# ========================================================================
# 3. Trace Matrices Precomputation
# ========================================================================

def comp_tracematrices_ramoddp(unif, d, e, maxpower, E, F, trace_func, basis, notnorm, anti=False):
    """
    Precompute trace matrices for ramified odd primes.
    
    INPUT:
        unif       : element; uniformizer of the local field F
        d, e       : integers; different valuation and ramification index
        maxpower   : integer; maximum depth parameter
        E, F       : algebraic fields; extension E and base F
        trace_func : function; trace map from E to F
        basis      : list; basis of E over F
        notnorm    : element; an element in F that is not a norm from E
        anti       : boolean; if True, computes for anti-Hermitian spaces

    OUTPUT:
        A list of precomputed trace matrices (or pairs of matrices for norm/non-norm variants).
    """
    L = []
    delta = E.gen() if anti else 1 
    
    for l in range(-d, e * maxpower + 1):
        if l % 2 == 0:
            M1 = matrix(E, 1, 1, [delta * unif**(l//2)])
            T1 = tracematrix(M1, E, F, trace_func, basis)
            M2 = matrix(E, 1, 1, [delta * unif**(l//2) * notnorm])
            T2 = tracematrix(M2, E, F, trace_func, basis)
            L.append([T1, T2])
        else:
            val = unif**((l+1)//2)
            M = matrix(E, 1, 1, [delta * val])
            T = tracematrix(M, E, F, trace_func, basis)
            L.append(block_matrix(F, 2, 2, [[0, T], [T.transpose(), 0]], subdivide=False))
    return L

def comp_tracematrices_ram2(unif, d, e, maxpower, E, F, trace_func, basis, anti=False):
    """
    Precompute trace matrices for ramified primes where p = 2.
    
    INPUT:
        unif       : element; uniformizer of the local field F
        d, e       : integers; different valuation and ramification index
        maxpower   : integer; maximum depth parameter
        E, F       : algebraic fields; extension E and base F
        trace_func : function; trace map from E to F
        basis      : list; basis of E over F
        anti       : boolean; if True, computes for anti-Hermitian spaces

    OUTPUT:
        A list of precomputed trace matrices tailored for the p=2 ramified case.
    """
    L = []
    delta = E.gen() if anti else 1
    for l in range(-d, e * maxpower + 1):
        if l % 2 == 1:
            val = unif**((l+1)//2)
            M = matrix(E, 1, 1, [delta * val])
            T = tracematrix(M, E, F, trace_func, basis)
            L.append([block_matrix(F, 2, 2, [[0, T], [T.transpose(), 0]], subdivide=False)])
        else:
            M1 = matrix(E, 1, 1, [delta * unif**(l//2)])
            T1 = tracematrix(M1, E, F, trace_func, basis)
            L.append([T1])
    return L

def comp_tracematrices_unr(unif, d, e, maxpower, E, F, trace_func, basis):
    """
    Precompute trace matrices for unramified cases.
    
    INPUT:
        unif       : element; uniformizer of the local field F
        d, e       : integers; different valuation (0) and ramification index (1)
        maxpower   : integer; maximum depth parameter
        E, F       : algebraic fields; extension E and base F
        trace_func : function; trace map from E to F
        basis      : list; basis of E over F

    OUTPUT:
        A list of trace matrices corresponding to unramified filtration steps.
    """
    matlist = []
    for l in range(-d, e * maxpower + 1):
        M = matrix(E, 1, 1, [unif**l])
        matlist.append(tracematrix(M, E, F, trace_func, basis))
    return matlist

# ========================================================================
# 4. Lattice Construction and Enumeration
# ========================================================================

def unique_hermi_unr(q, p, e, mult, E, OE, trace_func, basis):
    """
    Construct the unique maximal unramified Hermitian lattice.
    Used exclusively when maxpower == 0 (no p-power structure required).
    
    INPUT:
        q          : integer; residue field size of F
        p          : integer; prime number
        e          : integer; ramification index (always 1 here)
        mult       : integer; dimension/multiplicity
        E          : algebraic field; extension field
        OE         : integer ring; ring of integers of E
        trace_func : function; trace map from E to F
        basis      : list; basis of E over F

    OUTPUT:
        A list containing a single list `[isnorm, lattice_obj]` where `lattice_obj` 
        is an instantiated hermilattice_unr object.
    """
    isnorm = 0 
    
    A = identity_matrix(E, mult)
    gamma_mat = identity_matrix(mult * e)
    
    mass_term = hermi_unr_mass_local_term(q, mult, [mult])
    lattice_obj = hermilattice_unr(q, p, 0, mult * e, [mult], A, gamma_mat, mass_term)
    
    gram_F = tracematrix(A, E, OE.base_ring(), trace_func, basis)
    lattice_obj.Qpclass = quadQp_from_mat(gram_F, p)
    
    return [[isnorm, lattice_obj]]

def enum_hermi(q, p, unif, E, OE, residue_field, conj, val_func, trace_func, basis, 
              d, e, mult, maxpower, notnorm=None):
    """
    Enumerate all isomorphism classes of Hermitian lattices over a local field.
    Assumes the extension E/F is a field extension (inert or ramified).
    
    INPUT:
        q, p, unif : algebraic data; residue field size, prime characteristic, and uniformizer
        E, OE      : algebraic fields; extension field and its ring of integers
        residue_field, conj, val_func, trace_func, basis: extension operations/data
        d, e       : integers; different valuation and ramification index
        mult       : integer; total dimension
        maxpower   : integer; maximum depth of the p-power structure
        notnorm    : element (optional); an element in F not a norm from E (needed for odd ramified)

    OUTPUT:
        A list of `[invariant, lattice_obj]` containing all valid Hermitian lattices 
        and their associated invariant properties.
    """
    result = []
    
    # ---------------------------------------------------------------------
    # 1. Ramified cases (d > 0)
    # ---------------------------------------------------------------------
    if d > 0:
        if p == 2:
            formal_types = formal_types_hermi_ram2(d, e, mult, maxpower)
            p_t_matrices = comp_tracematrices_ram2(unif, d, e, maxpower, E, OE.base_ring(), trace_func, basis)
            
            for I in formal_types:
                # Append built ramified p=2 lattices here
                pass 
                
        else:
            formal_types = enum_hermi_ramoddp_types(d, e, mult, maxpower)
            if notnorm is None:
                raise ValueError("Ramified odd case requires 'notnorm' element.")
            
            p_t_matrices = comp_tracematrices_ramoddp(unif, d, e, maxpower, E, OE.base_ring(), 
                                                      trace_func, basis, notnorm)
            
            for I in formal_types:
                # Append built ramified odd p lattices here
                pass

    # ---------------------------------------------------------------------
    # 2. Unramified cases (d = 0)
    # ---------------------------------------------------------------------
    else:
        effective_mult = mult 
        
        if maxpower == 0:
            return unique_hermi_unr(q, p, e, effective_mult, E, OE, trace_func, basis)
        else:
            ftypes = partitions(effective_mult, maxpower)
            precomp_trace_matrices = comp_tracematrices_unr(unif, 0, e, maxpower, E, OE.base_ring(), trace_func, basis)
            
            for formal_type in ftypes:
                latt_data = build_hermilattice_unr(
                    q, p, unif, 0, 0, d, e, mult, formal_type, precomp_trace_matrices,
                    E, OE, residue_field, conj, val_func
                )
                
                gram_F = tracematrix(latt_data[1].gram_mat, E, OE.base_ring(), trace_func, basis)
                latt_data[1].Qpclass = quadQp_from_mat(gram_F, p)
                
                result.append(latt_data)
                
    return result
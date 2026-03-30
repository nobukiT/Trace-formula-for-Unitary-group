load("hermilatticeclass.sage") 
load("young.sage")
load("utils.sage")

# ========================================================================

def enum_partitions(n, maxpower):
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
        for K in enum_partitions(n - i, maxpower - 1):
            L.append([i] + K)
        i -= 1
    return L

def tracematrix(A, alg_data):
    """
    Compute the trace matrix of a Hermitian matrix A over the extension E/F.
    """
    Ei = alg_data['E']
    conj = alg_data['conj']
    
    # Ei.basis() might fail on some category objects, 
    # so we use integral_basis() which is more robust for number fields,
    # or access the basis via the absolute_field if available.
    try:
        basis_E = list(Ei.basis())
    except AttributeError:
        # Fallback for category objects or absolute fields
        basis_E = list(Ei.absolute_generator()**i for i in range(Ei.degree()))
    
    nrows, ncols = A.dimensions()
    deg = len(basis_E)
    matlist = []
    
    for a in range(nrows):
        for b in range(ncols):
            entries = [[QQ((A[a,b] * basis_E[i] * conj(basis_E[j])).trace()) 
                        for j in range(deg)] for i in range(deg)]
            matlist.append(matrix(QQ, deg, deg, entries))
            
    return block_matrix(nrows, ncols, matlist, subdivide=False)

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

def enum_hermi_ram2_types(d, e, mult, maxpower, mini=0):
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
            for J in enum_hermi_ram2_types(d, e, mult - 2 * i, maxpower, mini + 1):
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
            for J in enum_hermi_ram2_types(d, e, mult - i, maxpower, mini + 1):
                if i % 2 == 0:
                    L.append([[i, Mod(0, 2)]] + J)
                if i > 0:
                    L.append([[i, Mod(1, 2), Mod(0, 2)]] + J)
                    if J[0][0] == 0 and (len(J) == 1 or J[1][0] == 0 or J[1][1] == 0):
                        L.append([[i, Mod(1, 2), Mod(1, 2)]] + J)
    return L

# ========================================================================

def comp_tracematrices_ramoddp(maxpower, alg_data, anti=False):
    """
    Precompute trace matrices for ramified odd primes.
    
    INPUT:
        maxpower   : integer; maximum depth parameter
        alg_data
        anti       : boolean; if True, computes for anti-Hermitian spaces

    OUTPUT:
        A list of precomputed trace matrices (or pairs of matrices for norm/non-norm variants).
    """
    L = []
    delta = alg_data['E'].gen() if anti else 1 
    
    for l in range(-alg_data['d'], alg_data['e'] * maxpower + 1):
        if l % 2 == 0:
            M1 = matrix(alg_data['E'], 1, 1, [delta * alg_data['unif']**(l//2)])
            T1 = tracematrix(M1, alg_data)
            M2 = matrix(alg_data['E'], 1, 1, [delta * alg_data['unif']**(l//2) * alg_data['notnorm']])
            T2 = tracematrix(M2, alg_data)
            L.append([T1, T2])
        else:
            val = alg_data['unif']**((l+1)//2)
            M = matrix(alg_data['E'], 1, 1, [delta * val])
            T = tracematrix(M, alg_data)
            L.append(block_matrix(alg_data['F'], 2, 2, [[0, T], [T.transpose(), 0]], subdivide=False))
    return L

def comp_tracematrices_ram2(maxpower, alg_data, anti=False):
    """
    Precompute trace matrices for ramified primes where p = 2.
    
    INPUT:
        maxpower   : integer; maximum depth parameter
        alg_data
        anti       : boolean; if True, computes for anti-Hermitian spaces

    OUTPUT:
        A list of precomputed trace matrices tailored for the p=2 ramified case.
    """
    L = []
    delta = alg_data['E'].gen() if anti else 1
    for l in range(-alg_data['d'], alg_data['e'] * maxpower + 1):
        if l % 2 == 1:
            val = alg_data['unif']**((l+1)//2)
            M = matrix(alg_data['E'], 1, 1, [delta * val])
            T = tracematrix(M, alg_data)
            T_block = block_matrix(QQ, 2, 2, [[0, T], [T.transpose(), 0]], subdivide=False)
            L.append([T_block, T_block])
        else:
            M1 = matrix(alg_data['E'], 1, 1, [delta * alg_data['unif']**(l//2)])
            T1 = tracematrix(M1, alg_data)
            L.append([T1, T1, T1])
    return L

def comp_tracematrices_unr(maxpower, alg_data):
    """
    Precompute trace matrices for unramified cases.
    
    INPUT:
        maxpower   : integer; maximum depth parameter
        alg_data

    OUTPUT:
        A list of trace matrices corresponding to unramified filtration steps.
    """
    matlist = []
    for l in range(0, alg_data['e'] * maxpower + 1):  # Unramified so d=0
        M = matrix(alg_data['E'], 1, 1, [alg_data['unif']**l])
        matlist.append(tracematrix(M, alg_data))
    return matlist

# ========================================================================

def unique_hermi_unr(p, mult, alg_data):
    """
    Construct the unique maximal unramified Hermitian lattice.
    Used exclusively when maxpower == 0 (no p-power structure required).
    
    INPUT:
        p          : integer; prime number
        mult       : integer; dimension/multiplicity
        alg_data

    OUTPUT:
        A list containing a single list `[isnorm, lattice_obj]` where `lattice_obj` 
        is an instantiated hermilattice_unr object.
    """
    isnorm = 0 
    A = identity_matrix(alg_data['E'], mult)
    gram_F = tracematrix(A, alg_data)
    abs_dim = gram_F.nrows()
    
    gamma_mat = identity_matrix(abs_dim)
    
    mass_term = hermi_unr_mass_local_term(alg_data['q'], mult, [mult])
    lattice_obj = hermilattice_unr(alg_data['q'], p, 0, abs_dim, [mult], A, gamma_mat, mass_term)
    
    lattice_obj.Qpclass = quadQp_from_mat(gram_F, p)
    
    return [[isnorm, lattice_obj]]

def enum_hermi(p, m, D_K, mult, maxpower, alg_data):
    """
    Enumerate all isomorphism classes of Hermitian lattices over a local field.
    Assumes the extension E/F is a field extension (inert or ramified).

    INPUT:
        p          : integer; the prime number.
        m          : integer; the complete cyclotomic order.
        D_K        : integer; the discriminant of the global quadratic field.
        mult       : integer; total dimension (multiplicity).
        maxpower   : integer; maximum depth of the p-power structure.
        alg_data   : dict; dictionary containing local algebraic structures 
                     (q, unif, E, OE, residue_field, conj, val_func, trace_func, basis, d, e, notnorm).

    OUTPUT:
        list: A list of `[invariant, lattice_obj]` containing all valid Hermitian lattices.
    """
    k = Integer(m).valuation(p)
    m2 = m // (p**k)
    result = []
    
    if alg_data['d'] > 0:
        if p == 2:
            formal_types = enum_hermi_ram2_types(alg_data['d'], alg_data['e'], mult, maxpower)
            p_t_matrices = comp_tracematrices_ram2(maxpower, alg_data)
            for I in formal_types:
                latt = quadlattice_2_from_hermi_type(m, D_K, mult, I, p_t_matrices, alg_data)
                result.append(latt)
            return result
        else:
            formal_types = enum_hermi_ramoddp_types(alg_data['d'], alg_data['e'], mult, maxpower)
            if alg_data['notnorm'] is None:
                raise ValueError("Ramified odd case requires 'notnorm' element.")
            p_t_matrices = comp_tracematrices_ramoddp(maxpower, alg_data)
            for I in formal_types:
                latt = quadlattice_oddp_from_hermi_type(p, k, alg_data['d'], mult, I, p_t_matrices, alg_data)
                result.append(latt)
            return result
    else:
        effective_mult = mult 
        if maxpower == 0:
            return unique_hermi_unr(p, effective_mult, alg_data)
        else:
            ftypes = enum_partitions(effective_mult, maxpower)

            for formal_type in ftypes:
                latt_data = build_hermilattice_unr(
                    p, k, m2, mult, formal_type, alg_data
                )
                
                gram_F = tracematrix(latt_data[1].gram_mat, alg_data)
                latt_data[1].Qpclass = quadQp_from_mat(gram_F, p)
                
                result.append(latt_data)
        return result
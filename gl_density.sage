def gl_mass_local_term(q, I):
    """
    Calculate the local density (inverse mass) for a GL lattice over a local field F_v.
    INPUT:
        q : integer; the cardinality of the residue field of F_v (i.e., p^{f_v})
        I : list of integers; the partition type representing the lattice invariants
    OUTPUT:
        A rational number representing the local term in the mass formula.
    """
    n = sum(I)
    if n == 0:
        return 1
        
    mass_local_term = prod([(1 - q**(-l)) for l in range(1, int(n) + 1)])
    
    for i in range(len(I)):
        n_i = int(I[i])
        if n_i > 0:
            mass_local_term /= prod([(1 - q**(-l)) for l in range(1, n_i + 1)])
            for j in range(i + 1, len(I)):
                n_j = int(I[j])
                mass_local_term *= q**((j - i) * n_i * n_j)
                    
    return mass_local_term
load("young.sage")
load("local_density_ram2.sage")


def hermi_unr_mass_local_term(q, dim, I):
    """
    Calculate the local density of an unramified Hermitian lattice.
    INPUT:
        q       : integer; order of the residue field of the base field F_p
        dim     : integer; total dimension of the Hermitian space
        I       : list of integers; Jordan block sizes
    OUTPUT:
        A rational number representing the local density.
    """
    
    mass_local_term = prod([(1 - (-q)**(-l)) for l in range(1, int(dim) + 1)])
    for i in range(len(I)):
        n_i = int(I[i])
        if n_i > 0:
            mass_local_term /= prod([(1 - (-q)**(-l)) for l in range(1, n_i + 1)])
            for j in range(i + 1, len(I)):
                n_j = int(I[j])
                mass_local_term *= q**((j - i) * n_i * n_j)
            
    return mass_local_term

def hermi_ramoddp_mass_local_term(q, d, mult, I):
    """
    Calculate the local density of a ramified Hermitian lattice at an odd characteristic prime.

    INPUT:
        q    : integer; order of the residue field
        d    : integer; valuation of the different of the extension E/F
        mult : integer; total multiplicity (dimension) of the lattice
        I    : list of tuples; lattice invariants [(n_i, det_i), ...] where n_i is dimension
               and det_i is the determinant/type information.

    OUTPUT:
        A rational number representing the local density.
    """
    # Initial term based on the total dimension
    if mult % 2 == 0:
        mass_local_term = q**(mult // 2) * prod([1 - q**(-2*l) for l in range(1, 1 + int(mult // 2))])
    else:
        mass_local_term = prod([1 - q**(-2*l) for l in range(1, int((mult + 1) // 2))])
        
    q_exponent = 0
    for i in range(len(I)):
        n_i = int(I[i][0])
        
        if (i - d) % 2 == 0:
            if n_i > 0:
                if n_i % 2 == 1:
                    mass_local_term /= 2 * prod([1 - q**(-2*l) for l in range(1, int((n_i + 1) // 2))])
                else:
                    mass_local_term /= 2 * prod([1 - q**(-2*l) for l in range(1, int(n_i // 2))])
                    
                    cond1 = (I[i][1] % 2 == 0 and (q % 4 == 1 or n_i % 4 == 0))
                    cond2 = (I[i][1] % 2 == 1 and q % 4 == 3 and n_i % 4 == 2)
                    
                    if cond1 or cond2:
                        mass_local_term /= 1 - q**(-n_i // 2)
                    else:
                        mass_local_term /= 1 + q**(-n_i // 2)
                        
                for j in range(i + 1, len(I)):
                    n_j = int(I[j][0])  
                    q_exponent += (j - i) * n_i * n_j // 2
        else:
            if n_i > 0:
                mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int(n_i // 2))])
                q_exponent -= n_i // 2
                for j in range(i + 1, len(I)):
                    n_j = int(I[j][0])
                    q_exponent += (j - i) * n_i * n_j // 2
                    
    mass_local_term *= q**q_exponent
    return mass_local_term

def hermi_ram2_mass_local_term(Fi, Ei, embed_F_to_E, p_ideal, I):
    """
    Calculate the local mass term (inverse of Haar measure \nu_0) of a ramified 
    Hermitian lattice at a prime with residue characteristic 2.

    INPUT:
        Fi           : NumberField; Base totally real field
        Ei           : NumberField; CM field (quadratic extension of Fi)
        embed_F_to_E : Homomorphism; embedding from Fi to Ei
        p_ideal      : Prime ideal of Fi (assumed to be above 2)
        I            : list of tuples; lattice invariants 
                       (e.g., [n_i, type_i, ...] where n_i is dimension)

    OUTPUT:
        A rational number representing the local mass term (\nu_0^{-1}).
    """
    if Fi == QQ:
        p = p_ideal
        q = p
        e = 1
        p_in_E = Ei.ideal(p)
    else:
        q = p_ideal.norm()
        e = p_ideal.ramification_index()
        p_in_E = Ei.ideal([embed_F_to_E(g) for g in p_ideal.gens()])
    print(e)
    
    mult = sum(int(item[0]) if isinstance(item, (list, tuple)) else int(item) for item in I)
    
    ramified_primes = [(P, e_idx) for (P, e_idx) in p_in_E.factor() if e_idx > 1]
    
    if not ramified_primes:
        raise ValueError(f"Prime p is unramified in Ei. Cannot use ramified mass term formula. Factorization: {p_in_E.factor()}")
        
    P_ideal = ramified_primes[0][0]
    
    diff_Ei = Ei.different().valuation(P_ideal)
    if Fi == QQ:
        diff_Fi_lifted = 0
    else:
        diff_Fi_lifted = Ei.ideal([embed_F_to_E(g) for g in Fi.different().gens()]).valuation(P_ideal)
    d = diff_Ei - diff_Fi_lifted

    if e == 1:
        mass_local_term = prod([1 - q**(-2*l) for l in range(1, 1 + int(mult//2))])
        q_exponent = mult if mult % 2 == 0 else 0
        
        for i in range(len(I)):
            n_i = int(I[i][0])  
            
            if (i - d) % 2 == 1:
                actual_ni = 2 * n_i
                q_exponent -= n_i
                if actual_ni > 0:
                    if (i > 0 and I[i-1][1] % 2 == 1) or (i+1 < len(I) and I[i+1][1] % 2 == 1):
                        mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int(actual_ni//2))])
                    else:  
                        mass_local_term /= 2 * prod([1 - q**(-2*l) for l in range(1, int(actual_ni//2))])
                        if I[i][1] % 2 == 0:
                            mass_local_term /= 1 - q**(-actual_ni//2)
                        else:
                            mass_local_term /= 1 + q**(-actual_ni//2)
                    for j in range(i+1, len(I)):
                        n_j = int(I[j][0]) 
                        if (j - d) % 2 == 1:
                            q_exponent += (j-i) * actual_ni * n_j
                        else:
                            q_exponent += (j-i) * n_i * n_j
            elif n_i > 0:
                if I[i][1] % 2 == 0:
                    q_exponent -= n_i
                    mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int(n_i//2))])
                else:
                    if i+3 > len(I) or I[i+2][1] % 2 == 0:
                        q_exponent -= 1
                    mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int((n_i-1)//2))])
                for j in range(i+1, len(I)):
                    n_j = int(I[j][0])
                    if (j - d) % 2 == 1:
                        q_exponent += (j-i) * n_i * n_j
                    else:
                        q_exponent += (j-i) // 2 * n_i * n_j
                        
        mass_local_term *= q**q_exponent
        return mass_local_term

    else:
        print(f"Calculate brute-force local density for {I}")
        mass_local_term = prod([1 - q**(-2*l) for l in range(1, 1 + int(mult // 2))])
        mass_local_term *= q**(mult * sum(i * int(I[i][0]) for i in range(len(I))))

        if mult % 2 == 0:
            mass_local_term *= q**(mult * d / 2)
            
        mass_local_term /= calculate_local_density(Fi, Ei, embed_F_to_E, p_ideal, I, N=4)
        
        return mass_local_term
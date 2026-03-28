load("young.sage")

def hermi_unr_mass_local_term(q, dim, I, issplit):
    """
    Local density of an unramified Hermitian lattice (for general local fields)
    q: Order of the residue field of the base field F_p (q = p^f)
    dim: Dimension of the space
    I: List of Jordan block sizes
    issplit: Whether the extension is split (True/False)
    """
    epsilon = 1 if issplit else -1
    
    mass_local_term = prod([(1 - (epsilon*q)**(-l)) for l in range(1, int(dim) + 1)])
    for i in range(len(I)):
        mass_local_term /= prod([(1 - (epsilon*q)**(-l)) for l in range(1, int(I[i]) + 1)])
        for j in range(i+1, len(I)):
            mass_local_term *= q**((j-i)*I[i]*I[j])
            
    return mass_local_term

def hermi_ramoddp_mass_local_term(q, d, mult, I):
    """
    Local density of a branched Hermitian lattice at an odd characteristic
    q: Order of the residue field of the base field
    d: Value of the different
    mult: Multiplicity of the lattice
    I: Invariants of the lattice (Type and dimension)
    """
    if mult % 2 == 0:
        mass_local_term = q**(mult//2) * prod([1 - q**(-2*l) for l in range(1, 1 + int(mult//2))])
    else:
        mass_local_term = prod([1 - q**(-2*l) for l in range(1, int((mult+1)//2))])
        
    q_exponent = 0
    for i in range(len(I)):
        n_i = int(I[i][0])
        
        if (i - d) % 2 == 0:
            if n_i > 0:
                if n_i % 2 == 1:
                    mass_local_term /= 2 * prod([1 - q**(-2*l) for l in range(1, int((n_i+1)//2))])
                else:
                    mass_local_term /= 2 * prod([1 - q**(-2*l) for l in range(1, int(n_i//2))])
                    
                    # -1 が F_q で平方剰余かどうかの判定に q % 4 を使用
                    cond1 = (I[i][1] % 2 == 0 and (q % 4 == 1 or n_i % 4 == 0))
                    cond2 = (I[i][1] % 2 == 1 and q % 4 == 3 and n_i % 4 == 2)
                    
                    if cond1 or cond2:
                        mass_local_term /= 1 - q**(-n_i//2)
                    else:
                        mass_local_term /= 1 + q**(-n_i//2)
                        
                for j in range(i+1, len(I)):
                    n_j = int(I[j][0])  
                    q_exponent += (j-i) * n_i * n_j // 2
        else:
            if n_i > 0:
                mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int(n_i//2))])
                q_exponent -= n_i // 2
                for j in range(i+1, len(I)):
                    n_j = int(I[j][0])
                    q_exponent += (j-i) * n_i * n_j // 2
                    
    mass_local_term *= q**q_exponent
    return mass_local_term

def hermi_ram2_mass_local_term(q, d, e, mult, I):
    """
    Local density of the ramified Hermitian lattice at a prime point (residue characteristic 2)
    q: Order of the residue field of the base field
    d: Value of the different
    e: Absolute ramification index of the base field over Q₂
    mult: Multiplicity of the lattice
    I: Invariants of the lattice;
    """
    mass_local_term = prod([1 - q**(-2*l) for l in range(1, 1 + int(mult//2))])
    q_exponent = mult if mult % 2 == 0 else 0
    
    for i in range(len(I)):
        n_i = int(I[i][0])  
        
        if (i - d) % 2 == 1:
            actual_ni = 2 * n_i
            q_exponent -= n_i
            if actual_ni > 0:
                if (i > 0 and I[i-1][1] % 2 == 1) or (i+1 < len(I) and I[i+1][1] % 2 == 1):  # "bound"
                    mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int(actual_ni//2))])
                else:  # "free"
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
            if I[i][1] % 2 == 0:  # "type II"
                q_exponent -= n_i
                mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int(n_i//2))])
            else:  # "type I"
                # 注意: 一般の絶対分岐指数 e > 1 の場合、Choの公式では
                # この term "beta" (i+3 > len(I) などの条件) や Type I/II の干渉条件が
                # e に依存してより複雑にシフトします。
                # 現在の構造は e=1 の場合の Type I/II 構造を前提としています。
                if i+3 > len(I) or I[i+2][1] % 2 == 0:
                    q_exponent -= 1  # term "beta"
                mass_local_term /= prod([1 - q**(-2*l) for l in range(1, 1 + int((n_i-1)//2))])
            for j in range(i+1, len(I)):
                n_j = int(I[j][0])
                if (j - d) % 2 == 1:
                    q_exponent += (j-i) * n_i * n_j
                else:
                    q_exponent += (j-i) // 2 * n_i * n_j
                    
    mass_local_term *= q**q_exponent
    return mass_local_term
def gl_mass_local_term(p, k, m2, I):
    f = euler_phi(m2) if m2 > 2 else 1
    q = p**f
    n = sum(I)
    mass_local_term = prod([(1 - q**(-l)) for l in range(1, int(n) + 1)])
    for i in range(len(I)):
        n_i = int(I[i])
        if n_i > 0:
            mass_local_term /= prod([(1 - q**(-l)) for l in range(1, n_i + 1)])
            for j in range(i + 1, len(I)):
                n_j = int(I[j])
                mass_local_term *= q**((j - i) * n_i * n_j)
                    
    return mass_local_term
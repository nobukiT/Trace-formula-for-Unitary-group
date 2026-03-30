def calculate_local_density(Fi, Ei, embed_F_to_E, p_ideal, I_local, N=2):
    """
    Fi : base field (here Q(sqrt(2)))
    Ei : absolute CM field containing Fi (here Q(zeta_8))
    embed_F_to_E : embedding Fi -> Ei
    p_ideal : prime ideal of Fi above 2
    I_local : Jordan data [(n_i, type_i), ...]
    N : modulus level

    Empirical brute-force local density over O_E / P^(2N).
    """
    if Fi == QQ:
        p = p_ideal
        q = p
        p_in_E = Ei.ideal(p)
        unif_F = QQ(p)
        unif_E = Ei(p)
    else:
        q = p_ideal.norm()
        p_in_E = Ei.ideal([embed_F_to_E(g) for g in p_ideal.gens()])
        gens_F = p_ideal.gens_reduced()
        unif_F = None
        p_sq = p_ideal**2
        for g in gens_F:
            if g != 0 and g not in p_sq:
                unif_F = Fi(g)
                break
        if unif_F is None:
            raise ValueError("Could not find unif_F.")
        unif_E = embed_F_to_E(unif_F)

    factored_primes = list(p_in_E.factor())
    ramified_primes = [(P, e) for (P, e) in factored_primes if e > 1]
    P_ideal = ramified_primes[0][0]

    n_dim = sum(n for n, t in I_local)
    if n_dim == 0:
        return 1

    max_power = N

    conj_func = None
    for aut in Ei.automorphisms():
        if aut(Ei.gen()) != Ei.gen():
            if Fi == QQ:
                conj_func = aut
                break
            else:
                if all(aut(embed_F_to_E(g)) == embed_F_to_E(g) for g in Fi.gens()):
                    conj_func = aut
                    break
    if conj_func is None:
        conj_func = Ei.complex_conjugation()

    gens_P = P_ideal.gens_reduced()
    unif_P = None
    P_sq = P_ideal**2
    for g in gens_P:
        if g != 0 and g not in P_sq:
            unif_P = Ei(g)
            break
    if unif_P is None:
        raise ValueError("Could not find unif_P.")

    OE = Ei.ring_of_integers()
    OE_mod = OE.quotient(P_ideal**max_power, names='ubar')

    def conj_mod(x_mod):
        x_lift = Ei(x_mod.lift())
        x_conj = conj_func(x_lift)
        return OE_mod(x_conj)

    blocks = []
    for i, (n_i, type_i) in enumerate(I_local):
        if n_i == 0:
            continue

        scale_val = unif_E**i
        scale_mod = OE_mod(scale_val)

        if type_i % 2 == 0:  # Type II
            if n_i % 2 != 0:
                raise ValueError("Type II blocks must have even dimension.")
            H = matrix(OE_mod, 2, 2, [0, scale_mod, conj_mod(scale_mod), 0])
            blocks.extend([H] * (n_i // 2))
        else:  # Type I
            D = matrix(OE_mod, 1, 1, [scale_mod])
            blocks.extend([D] * n_i)

    B_mat = block_diagonal_matrix(blocks)
    
    k, k_map = Ei.residue_field(P_ideal)
    k_reps = [Ei(k_map.section()(a)) for a in k]
    
    ring_elements = []
    for coeffs in product(k_reps, repeat=max_power):
        val = sum(c * (unif_P**j) for j, c in enumerate(coeffs))
        ring_elements.append(OE_mod(val))

    valid_count = 0
    for matrix_entries in product(ring_elements, repeat=n_dim*n_dim):
        X = matrix(OE_mod, n_dim, n_dim, matrix_entries)
        X_star = X.apply_map(conj_mod).transpose()

        if X_star * B_mat * X == B_mat:
            valid_count += 1

    dim_G = n_dim**2
    empirical_density = (q**(-N * dim_G)) * valid_count

    return empirical_density
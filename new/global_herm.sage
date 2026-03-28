load("Lfunc.sage")
load("orb_int_common.sage")


# ========================================================================
def inverse_root_polynomial(f):
    """
    Return x^d * f(1/x), monic.
    """
    R = f.parent()
    x = R.gen()
    d = f.degree()
    return (x**d * f(x**(-1))).expand().monic()


def factor_dict(poly):
    return {g.monic(): e for g, e in poly.factor()}


def is_elliptic(m, K):
    """
    Return True iff Phi_m(x), viewed over K[x], has NO GL-type obstruction.
    That is: every irreducible factor f of Phi_m over K is paired with
    inverse_root_polynomial(f) with the same multiplicity.
    Since Phi_m is squarefree, this means:
      for every irreducible factor f, inverse_root_polynomial(f) is also
      a factor of Phi_m over K.
    This is the correct criterion for "Phi_m contributes no GL-type factor".
    """
    R.<x> = PolynomialRing(K)
    phi = R(cyclotomic_polynomial(m))
    fac = factor_dict(phi)

    for f in fac:
        finv = inverse_root_polynomial(f)
        if finv not in fac:
            return False
    return True

def append_poss_local_invs_unitary(poss_local_invs, p, sub_cc, database, D_K, verbose=False):
    """
    For a fixed prime p and a grouped local decomposition sub_cc,
    append all possible p-local invariants to poss_local_invs.

    INPUT:
        poss_local_invs : list to which the p-local possibilities are appended
        p               : prime
        sub_cc          : output of local_p_subclasses_unitary(cc, p)
        database        : cache / lookup table
        D_K             : discriminant of the quadratic field K
    OUTPUT:
        poss_local_invs.append(p_invs)
        where p_invs is a list of:
            [HW_list_ordered, orb_int]
    """
    p_invs = []
    block_choices = []

    for block in sub_cc:
        m_red, dims_tup, indices = block
        db_key = (p, m_red, tuple(dims_tup), D_K)

        if db_key not in database:
            if verbose:
                print(f"Warning: cache miss for {db_key}. This conjugacy class drops to Mass 0.")
            return False

        orb_dict = database[db_key]
        current_block_choices = []

        for invs_tup, orb_int in orb_dict.items():
            HW_list = [invs_tup[i] for i in range(len(indices))]
            current_block_choices.append([HW_list, orb_int])

        if len(current_block_choices) == 0:
            return False

        block_choices.append(current_block_choices)

    total_length = sum(len(block[2]) for block in sub_cc)

    for choice_tuple in itertools.product(*block_choices):
        HW_list_ordered = [None] * total_length
        orb_int = 1

        for block_idx, block in enumerate(sub_cc):
            _, _, indices = block
            sub_HW_list, sub_orb_int = choice_tuple[block_idx]

            orb_int *= sub_orb_int

            for j, ind in enumerate(indices):
                HW_list_ordered[ind] = sub_HW_list[j]

        p_invs.append([HW_list_ordered, orb_int])

    poss_local_invs.append(p_invs)
    return True

# ========================================================================
# 3. Calculation of mass
# ========================================================================
def mass_global_term_unitary(cc, D_K, prec=100, verbose=False):
    """
    Compute the global blockwise term for a unitary centralizer attached to
        cc = [[m1,n1],[m2,n2],...].
    INPUT:
        cc      -- list [[m,n],...]
        D_K     -- quadratic field parameter for E = Q(sqrt(D_K))
        prec    -- precision for L-values
        verbose -- print debugging info
    OUTPUT:
        product over blocks of the prescribed global factors
    """
    CC = ComplexField(prec)
    res = CC(1)

    for m, n in cc:
        data = block_fields_unitary(m, D_K, verbose)
        Ei = data["Ei"]
        Fi = data["Fi"]
        deg_Fi = data["deg_Fi"]
        
        if Fi is None:
            raise ValueError(f"Could not determine Fi for Phi_{m}-block.")
        m_red = m // 2 if m % 4 == 2 else m

        # ---- Tamagawa number ----
        res *= 2

        # ---- L-factor ----
        res *= unitary_block_L_product(n, Ei, Fi, prec, verbose)

        if verbose:
            print("=" * 80)
            print(f"Block Phi_{m}^{n}")
            print(f"Ei = {Ei}")
            print(f"Fi = {Fi}")
            print(f"block_size = {n}")
            print(f"partial result = {res}")
            print("=" * 80)

    return res

def mass_unitary(n, primes_list, negdim_goal, cc, database, D_K, prec=100, verbose=False):
    poss_local_invs = []

    # ---- finite places ----
    for p in primes_list:
        sub_cc = tot_ram_sub_conj_classes_unitary(cc, p)
        if not append_poss_local_invs_unitary(poss_local_invs, p, sub_cc, database, D_K):
            return 0

    mass = 0
    num_primes = len(primes_list)

    # Precompute field data for each Phi_m-block
    block_info = []
    for m, mult in cc:
        data = block_fields_unitary(m, D_K, False)
        block_info.append({
            "m": m,
            "mult": mult,
            "Ei": data["Ei"],
            "Fi": data["Fi"],
            "deg_Fi": data["deg_Fi"]
        })

    # ---- sum over finite local invariants ----
    for finite_invs in itertools.product(*poss_local_invs):
        orb_int = prod([finite_invs[i][1] for i in range(num_primes)])

        negdim_list = []

        for j, B in enumerate(block_info):
            m = B["m"]
            mult = B["mult"]
            num_real_places = B["deg_Fi"]

            HW_sum = sum([int(finite_invs[i][0][j]) for i in range(num_primes)]) % 2

            poss_negdims = []
            for archi_negdims in itertools.product(range(mult + 1), repeat=num_real_places):
                if sum(archi_negdims) % 2 == HW_sum:
                    arch_vol = prod([
                        (-1)**(negdim * (mult - negdim)) * binomial(mult, negdim) / (2**mult)
                        for negdim in archi_negdims
                    ])
                    total_negdim = sum(archi_negdims)
                    poss_negdims.append([total_negdim, arch_vol])

            negdim_list.append(poss_negdims)

        # ---- combine archimedean choices across blocks ----
        for archi_negdims in itertools.product(*negdim_list):
            if sum([x[0] for x in archi_negdims]) == negdim_goal:
                mass += orb_int * prod([x[1] for x in archi_negdims])

    global_factor = mass_global_term_unitary(cc, D_K, prec=prec, verbose=verbose)
    return mass * global_factor
# ========================================================================
# 4. main function
# ========================================================================
def mass_list_unitary(n, primes_list, database, D_K, imposed_negdim=-1):
    classes = conjclasses_unitary(n, primes_list)
    negdim_goal = n if imposed_negdim == -1 else imposed_negdim
    
    mass_list = []
    print(f"Starting calculation for U({n},{n}) with D_K = {D_K}")
    
    for conj_class in classes:
        cc_mass = mass_unitary(n, primes_list, negdim_goal, conj_class, database, D_K)
        
        if cc_mass != 0:
            mass_list.append([conj_class, cc_mass])
            print(f"{conj_class} -> Mass: {cc_mass}")
            
    print(f"Found {len(mass_list)} conjugacy classes with non-zero mass.")
    return mass_list
load("Lfunc.sage")
load("orb_int_common.sage")
load("utils.sage")

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
# Calculation of mass
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
# main function
# ========================================================================
def mass_list_unitary(n, primes_list, database, D_K, imposed_negdim=-1, prec=100, verbose=False):
    """
    Enumerate stable conjugacy classes of finite-order elliptic semisimple elements
    in U(p,q) (with p+q = n in the current normalization), and compute the mass
    attached to each class.

    INPUT:
        n               -- integer; we are working with total E-dimension n
        primes_list     -- list of primes used for local ramification / database lookup
        database        -- precomputed local invariant database
        D_K             -- discriminant parameter for E = Q(sqrt(D_K))
        imposed_negdim  -- target q for the original U(p,q)
                           if -1, default to q = n/2 (i.e. U(n/2,n/2))
        prec            -- precision for global L-values
        verbose         -- print detailed debugging information

    OUTPUT:
        list of pairs
            [conj_class, mass]
        where
            conj_class = [[m1,a1],[m2,a2],...]
        means characteristic polynomial
            prod_i Phi_{m_i}(x)^{a_i}
        and 'mass' is the total contribution from that stable class.
    """

    # ------------------------------------------------------------
    # Enumerate all admissible elliptic characteristic polynomials
    # of total Q-dimension n.
    # ------------------------------------------------------------
    classes = conjclasses_unitary_elliptic(n, primes_list)

    # ------------------------------------------------------------
    # Choose the target negative dimension q of the original U(p,q).
    # Default: q = n/2.
    # ------------------------------------------------------------
    negdim_goal = n if imposed_negdim == -1 else imposed_negdim

    if negdim_goal < 0 or negdim_goal > n:
        raise ValueError(f"imposed_negdim must satisfy 0 <= q <= {n}, got {negdim_goal}")

    mass_list = []

    print("=" * 80)
    print(f"Starting mass computation for unitary group with total size {n}")
    print(f"Quadratic field parameter D_K = {D_K}")
    print(f"Target signature parameter q = {negdim_goal}")
    print(f"Number of candidate stable classes = {len(classes)}")
    print("=" * 80)

    # ------------------------------------------------------------
    # 3. Loop over all admissible stable conjugacy classes.
    #    Each conj_class is a list [[m1,a1],[m2,a2],...]
    #    corresponding to charpoly prod Phi_m(x)^a.
    # ------------------------------------------------------------
    for conj_class in classes:
        if verbose:
            print("-" * 80)
            print(f"Processing class: {conj_class}")

        cc_mass = mass_unitary(
            n=n,
            primes_list=primes_list,
            negdim_goal=negdim_goal,
            cc=conj_class,
            database=database,
            D_K=D_K,
            prec=prec,
            verbose=verbose
        )

        # Keep only nonzero-mass classes
        if cc_mass != 0:
            mass_list.append([conj_class, cc_mass])
            print(f"{conj_class} -> Mass: {cc_mass}")
        elif verbose:
            print(f"{conj_class} -> Mass 0 (discarded)")

    print("=" * 80)
    print(f"Found {len(mass_list)} stable conjugacy classes with non-zero mass.")
    print("=" * 80)

    return mass_list
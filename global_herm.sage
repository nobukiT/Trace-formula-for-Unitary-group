load("Lfunc.sage")
load("utils.sage")
import itertools

def append_poss_local_invs_unitary(poss_local_invs, p, sub_cc, orb_int_unitary_db, D_K, total_num_blocks):
    """
    Fetches all valid local invariant patterns and their densities from the database
    for a given prime p and its local conjugacy class decomposition.

    INPUT:
        poss_local_invs    : list; a list to which the local possibilities will be appended.
        p                  : integer; the prime number.
        sub_cc             : list; the output of tot_ram_sub_conj_classes_unitary(cc, p).
        orb_int_unitary_db : dict; the precomputed local orbital integral database.
        D_K                : integer; the discriminant of the quadratic field E = Q(sqrt(D_K)).

    OUTPUT:
        Returns True if all required local data exist in the DB and are appended.
        Otherwise returns False (indicating mass 0 for this global class).
        Format appended to poss_local_invs:
            [ ((inv1), (inv2), ...), density ]
            where (inv_i) is the local invariant (e.g., (rank, disc)) for each Jordan block.
    """
    p_local_options = []
    
    for cyclo_order, jordan_blocks, orig_indices in sub_cc:
        key = (p, cyclo_order, tuple(jordan_blocks), D_K)
        if key not in orb_int_unitary_db:
            return False
        
        local_map = orb_int_unitary_db[key]
        # Store (mapped_pattern, density)
        # mapped_pattern is a dict {global_idx: local_inv}
        options = []
        for pattern, density in local_map.items():
            mapped_pattern = {orig_indices[k]: pattern[k] for k in range(len(orig_indices))}
            options.append((mapped_pattern, density))
        p_local_options.append(options)

    for combo in itertools.product(*p_local_options):
        combined_pattern_dict = {}
        combined_density = QQ(1)
        for p_map, density in combo:
            combined_pattern_dict.update(p_map)
            combined_density *= density
        
        # Convert dict to a list ordered by global block index
        ordered_pattern = [combined_pattern_dict[i] for i in range(total_num_blocks)]
        poss_local_invs.append([ordered_pattern, combined_density])
    
    return True

def mass_global_term_unitary(cc, D_K, prec=100, verbose=False):
    """
    Computes the product of global factors (Tamagawa number, L-values) for a unitary centralizer.

    INPUT:
        cc      : list; the stable conjugacy class [[m1, n1], [m2, n2], ...].
        D_K     : integer; discriminant of the quadratic field E = Q(sqrt(D_K)).
        prec    : integer; precision for L-value calculation.
        verbose : boolean; print debug information.

    OUTPUT:
        The product of global constants and L-factors over all blocks.
    """
    CC = ComplexField(prec)
    res = CC(1)

    for m, n in cc:
        data = block_fields_unitary(m, D_K, verbose)
        Ei = data["Ei"]
        Fi = data["Fi"]
        
        if Fi is None:
            raise ValueError(f"Could not determine Fi for Phi_{m}-block.")

        # Tamagawa number factor
        res *= 2

        # L-factor calculation
        res *= unitary_block_L_product(n, Ei, Fi, verbose)

    return res

def mass_unitary(n, primes_list, negdim_goal, cc, database, D_K, prec=100, verbose=False):
    """
    Computes the total mass for a single stable conjugacy class.

    INPUT:
        n           : integer; total dimension of the unitary space.
        primes_list : list; primes to check for local consistency.
        negdim_goal : integer; target signature q in U(p,q) where p+q=n.
        cc          : list; the conjugacy class [[m, mult], ...].
        database    : dict; the local orbital integral database.
        D_K         : integer; discriminant of E.
        prec        : integer; precision for complex numbers.

    OUTPUT:
        A numerical value (mass) representing the weighted number of lattices 
        fixed by an element in this stable class.
    """
    total_num_blocks = len(cc)
    poss_local_invs_per_prime = []

    for p in primes_list:
        p_invs = []
        sub_cc = tot_ram_sub_conj_classes_unitary(cc, p)
        if not append_poss_local_invs_unitary(p_invs, p, sub_cc, database, D_K, total_num_blocks):
            return 0
        poss_local_invs_per_prime.append(p_invs)

    total_mass = 0
    num_primes = len(primes_list)
    block_info = []
    for m, mult in cc:
        data = block_fields_unitary(m, D_K, False)
        block_info.append({"mult": mult, "deg_Fi": data["deg_Fi"]})

    for finite_invs in itertools.product(*poss_local_invs_per_prime):
        orb_int = prod([x[1] for x in finite_invs])
        if orb_int == 0: continue

        negdim_options_per_block = []
        for j, B in enumerate(block_info):
            mult = B["mult"]
            num_real_places = B["deg_Fi"]
            HW_sum = sum([int(finite_invs[i][0][j]) for i in range(num_primes)]) % 2
            poss_negdims_for_block = []

            for archi_negdims in itertools.product(range(mult + 1), repeat=num_real_places):
                if sum(archi_negdims) % 2 == HW_sum:
                    arch_vol = prod([
                        (-1)**(negdim * (mult - negdim)) * binomial(mult, negdim) / (2**mult)
                        for negdim in archi_negdims
                    ])
                    poss_negdims_for_block.append([sum(archi_negdims), arch_vol])
            
            if not poss_negdims_for_block: return 0
            negdim_options_per_block.append(poss_negdims_for_block)

        for block_choices in itertools.product(*negdim_options_per_block):
            if sum([x[0] for x in block_choices]) == negdim_goal:
                total_mass += orb_int * prod([x[1] for x in block_choices])

    global_factor = mass_global_term_unitary(cc, D_K, prec=prec, verbose=verbose)
    return total_mass * global_factor

def mass_list_unitary(n, primes_list, database, D_K, imposed_negdim=-1, prec=100, verbose=False):
    """
    Main entry point for calculating the mass list of a unitary group.

    INPUT:
        n              : integer; the total dimension.
        primes_list    : list; primes considered.
        database       : dict; local database.
        D_K            : integer; discriminant.
        imposed_negdim : integer; the signature 'q' (default n).
        prec           : integer; computation precision.

    OUTPUT:
        list; [[conj_class, mass], ...] for all classes with non-zero mass.
    """
    classes = conjclasses_unitary(n, primes_list, D_K, verbose)
    negdim_goal = n//2 if imposed_negdim == -1 else imposed_negdim

    if negdim_goal < 0 or negdim_goal > n:
        raise ValueError(f"imposed_negdim must be 0 <= q <= {n}")

    mass_list = []
    for conj_class in classes:
        cc_mass = exactify_to_QQ(mass_unitary(n, primes_list, negdim_goal, conj_class, database, D_K, prec, verbose))
        if cc_mass != 0:
            mass_list.append([conj_class, cc_mass])
            print(f"{conj_class} -> Mass: {cc_mass}")

    return mass_list
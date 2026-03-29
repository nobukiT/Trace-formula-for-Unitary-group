load("orb_int_common.sage")
load("gl.sage")
load("utils.sage")
load("hermi.sage")
load("jordan.sage")
load("Lfunc.sage")
load("quadlatticeclass.sage")
load("hermilatticeclass.sage") 
load("gllatticeclass.sage") 

if 'first_run' not in locals():
    first_run = True

if first_run:
    orb_int_unitary_db = {}
    hermi_latt_db = {}  
    gl_latt_db = {}
    first_run = False

# ========================================================================
# 格子の列挙 (一般化された関数を使用)
# ========================================================================

def enum_with_inv_hermi(p, m, mult, maxpower, D_K, alg_data):
    """一般化された引数を用いてエルミート格子を列挙する。"""
    # 前述の一般化された enum_hermi を呼び出し
    Lattices = enum_hermi(
        alg_data['q'], p, alg_data['unif'], alg_data['E'], alg_data['OE'], 
        alg_data['residue_field'], alg_data['conj'], alg_data['val_func'], 
        alg_data['trace_func'], alg_data['basis'], alg_data['d'], alg_data['e'], 
        mult, maxpower
    )
    
    result = {}
    for B in Lattices:
        # B = [isnorm, lattice_obj]
        lattice = B[1]
        dims = tuple(lattice.type_list)
        if dims not in result:
            result[dims] = []
        result[dims].append(lattice)
    return result

# ========================================================================
# 3. 軌道積分エンジン (正規化の重複を排除)
# ========================================================================

def orb_int_gl(args, Fi=None):
    """
    This engine bridge matches the call: 
    for _, val in orb_int_gl(single_block_args, Fi):
    
    args: list of (p, exponent, multiplicity, max_power_param)
    """
    results = []
    for p, exponent, multiplicity, max_power_param in args:
        # Determine k (p-primary) and m2 (tame) from exponent
        # Note: cyclo_order from the outer loop is passed via logic
        # Here we assume the setup mirrors the logic in the unitary DB update
        
        # In the context of the unitary DB, we need the tame part m2.
        # This should ideally be passed in or calculated from a global state.
        # For this bridge, we assume m2 is handled or fixed.
        m2 = 1 # Placeholder: should be derived from cyclo_order
        
        # Calculate e (ramification index of cyclotomic field at p)
        if exponent == 0:
            e = 1
        else:
            e = p**(exponent - 1) * (p - 1)
            
        max_slots = e * max_power_param
        partitions = enum_gl_partitions(multiplicity, max_slots)
        
        for I in partitions:
            val = gl_mass_local_term(p, exponent, m2, I)
            # Yielding a pair (partition, value) as expected by the loop
            results.append((tuple(I), val))
            
    return results

def orb_int_unitary(args_list, hermi_latt_db, D_K):
    for args in args_list:
        p, cyclo_order, jordan_blocks, partition_tup = args
        total_orbital_integral = 1
        max_power_param = len(jordan_blocks) - 1
        
        for i, (exponent, multiplicity) in enumerate(jordan_blocks):
            complete_order = cyclo_order * (p**exponent)
            lattice_query_key = (p, complete_order, multiplicity, max_power_param, D_K)
            target_partition = partition_tup[i]
            
            matched_lattices = hermi_latt_db[lattice_query_key][target_partition]
            # ここでも既に正規化済みの mass_local_term をそのまま掛ける
            total_orbital_integral *= matched_lattices[0].mass_local_term
            
        yield args, total_orbital_integral

# ========================================================================
# 4. データベースの更新処理 (メインループ)
# ========================================================================

def update_orb_int_unitary_db(n, primes_list, D_K, pretend=False, verbose=True):
    """
    Populate orb_int_unitary_db with all local orbital integral data needed
    for elliptic finite-order stable conjugacy classes in the unitary setting.

    INPUT:
        n           : integer; total size (default global setting U(n/2,n/2))
        primes_list : list of rational primes at which local data may be needed
        D_K         : integer; discriminant parameter for E = Q(sqrt(D_K))
        pretend     : boolean; if True, only print what would be computed, do not store
        verbose     : boolean; print progress

    OUTPUT:
        None (updates the global dictionary orb_int_unitary_db in place)
    """
    
    conj_classes = conjclasses_unitary(n, primes_list)
    required_params = {
        (p, cyclo_order, tuple(jordan_blocks))
        for p in primes_list
        for cc in conj_classes
        for cyclo_order, jordan_blocks, _ in tot_ram_sub_conj_classes_unitary(cc, p)
    }

    new_params = [
        param for param in required_params
        if (param[0], param[1], param[2], D_K) not in orb_int_unitary_db
    ]
    num_tasks = len(new_params)
    
    if verbose:
        print("=" * 80)
        print(f"Updating unitary orbital-integral DB for D_K = {D_K}")
        print(f"Number of required local tasks: {num_tasks}")
        print("=" * 80)

    for idx, (p, cyclo_order, jordan_blocks) in enumerate(sorted(new_params), 1):
        db_storage_key = (p, cyclo_order, tuple(jordan_blocks), D_K)

        if verbose:
            print("-" * 80)
            print(f"[{idx}/{num_tasks}] Computing local data for key:")
            print(f"    p = {p}")
            print(f"    cyclo_order = {cyclo_order}")
            print(f"    jordan_blocks = {jordan_blocks}")
            print(f"    D_K = {D_K}")
            
        if pretend:
            continue

        block_data = block_fields_unitary(cyclo_order, D_K, False)
        Ei = block_data["Ei"]
        Fi = block_data["Fi"]
        
        # Checking splitting behavior over the first prime above p for classification
        P_Fi = Fi.primes_above(p)[0]
        is_split = (len(Ei.primes_above(P_Fi)) == 2)
        
        local_results_map = {}
        max_power_param = len(jordan_blocks) - 1

        # ============================================================
        # Local place splits completely in Ei/Fi (GL case)
        # ============================================================
        if is_split:
            primes_v_above_p = Fi.primes_above(p)
            total_density_at_p = 1
            
            for v in primes_v_above_p:
                ev = v.ramification_index()
                fv = v.residue_class_degree()
                qv = p ** fv
                
                density_sums_per_block = []

                for exponent, multiplicity in jordan_blocks:
                    complete_order = cyclo_order * (p ** exponent)
                    gl_key = (v, complete_order, multiplicity, max_power_param)

                    if gl_key not in gl_latt_db:
                        # enum_split completely handles the absolute_e and relative_e math internally
                        # and returns a list of [disc_dummy, lattice_obj]
                        gl_latt_db[gl_key] = enum_split(
                            p, complete_order, multiplicity, max_power_param, ev, qv
                        )

                    lattices_in_block = gl_latt_db[gl_key]
                    block_densities = []

                    # Iterate over the generated lattice objects directly
                    for disc_dummy, latt_obj in lattices_in_block:
                        # The mass was already calculated and attached to the object during build_gl_lattice
                        val = latt_obj.mass_local_term
                        block_densities.append(val)

                    density_sums_per_block.append(sum(block_densities))

                density_at_v = prod(density_sums_per_block)
                total_density_at_p *= density_at_v

            if total_density_at_p != 0:
                dummy_inv_key = tuple(Mod(0, 2) for _ in jordan_blocks)
                local_results_map[dummy_inv_key] = total_density_at_p

        # ============================================================
        # Local place is inert or ramified in Ei/Fi (Hermitian case)
        # ============================================================
        else:
            alg_data = setup_local_algebraic_data(p, cyclo_order, D_K)
            partition_options_per_block = []

            for exponent, multiplicity in jordan_blocks:
                complete_order = cyclo_order * (p ** exponent)
                hermi_key = (p, complete_order, multiplicity, max_power_param, D_K)

                if hermi_key not in hermi_latt_db:
                    hermi_latt_db[hermi_key] = enum_with_inv_hermi(
                        p, complete_order, multiplicity, max_power_param, D_K, alg_data
                    )

                partition_options_per_block.append(list(hermi_latt_db[hermi_key].keys()))

            for partition_combination in itertools.product(*partition_options_per_block):
                calculation_args = [(p, cyclo_order, jordan_blocks, partition_combination)]

                for result_args, density_value in orb_int_unitary(
                    calculation_args, hermi_latt_db, D_K
                ):
                    if density_value != 0:
                        inv_key = result_args[3]
                        local_results_map[inv_key] = (
                            local_results_map.get(inv_key, 0) + density_value
                        )

        orb_int_unitary_db[db_storage_key] = local_results_map

        if verbose:
            print(f"Stored {len(local_results_map)} local invariant patterns.")

    if verbose:
        print("=" * 80)
        print(f"Database update complete for D_K = {D_K}.")
        print("=" * 80)
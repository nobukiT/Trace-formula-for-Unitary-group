load("orb_int_common.sage")
load("gl.sage")
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
# 1. 局所拡張の代数的パラメータのセットアップ
# ========================================================================

def setup_local_algebraic_data(p, m, D_K):
    """
    指定された p, m, D_K に対して、計算に必要なすべての代数オブジェクトを構築する。
    """
    # 拡大体 E = Q(zeta_m) の構成
    E_global = CyclotomicField(m, 'z')
    ideal_p = E_global.primes_above(p)[0]
    
    # 局所的なパラメータ
    q = ideal_p.norm()
    e_E = ideal_p.absolute_ramification_index()
    residue_field = E_global.residue_field(ideal_p)
    
    # 対合 (conj) とトレース関数の定義
    z = E_global.gen()
    conj = E_global.hom([z**(-1)])
    trace_func = lambda x: E_global(x + conj(x))
    val_func = lambda x: E_global(x).valuation(ideal_p) // e_E
    
    # 基礎体上の基底 (簡易的に [1] としているが、次数に合わせて拡張が必要)
    basis = [E_global(1)] 
    unif = E_global(p) # 基礎体が Q の場合は p
    
    # D_K によるさらなる拡大 (L/E) のチェック
    R = PolynomialRing(E_global, 'x')
    poly = R.gen()**2 - D_K
    
    if poly.is_irreducible():
        L_global = E_global.extension(poly, 'a')
        ideal_L = L_global.primes_above(ideal_p)[0]
        e_L = ideal_L.absolute_ramification_index()
        d_L = L_global.absolute_different().valuation(ideal_L)
    else:
        e_L = e_E
        d_L = E_global.different().valuation(ideal_p)

    return {
        'q': q, 'unif': unif, 'E': E_global, 'OE': E_global.ring_of_integers(),
        'residue_field': residue_field, 'conj': conj, 'val_func': val_func,
        'trace_func': trace_func, 'basis': basis, 'd': d_L, 'e': e_L, 'e_E': e_E
    }

# ========================================================================
# 2. 格子の列挙 (一般化された関数を使用)
# ========================================================================

def enum_with_inv_hermi(p, m, mult, maxpower, D_K, alg_data):
    """一般化された引数を用いてエルミート格子を列挙する。"""
    # 前述の一般化された enumhermi を呼び出し
    Lattices = enumhermi(
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

def orb_int_gl_engine(args_list):
    for args in args_list:
        p, m, dims_tup, invs_tup = args
        orbital_integral = 1
        m2 = m
        while m2 % p == 0:
            m2 //= p
        for i, (k_i, mult_i) in enumerate(dims_tup):
            # 既に L 因子で割られている gl_mass_local_term をそのまま使用
            raw_val = gl_mass_local_term(p, k_i, m2, invs_tup[i])
            orbital_integral *= raw_val
        yield args, orbital_integral

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

def update_orb_int_unitary_db(n, primes_list, D_K, pretend=False):
    unitary_classes = conjclasses_unitary(n, primes_list)
    required_params = {
        (p, cyclo_order, tuple(jordan_blocks))
        for p in primes_list
        for u_class in unitary_classes
        for cyclo_order, jordan_blocks, _ in tot_ram_sub_conj_classes_unitary(u_class, p)
    }
    
    new_params = [p for p in required_params if (p[0], p[1], p[2], D_K) not in orb_int_unitary_db]
    num_tasks = len(new_params)

    for idx, (p, cyclo_order, jordan_blocks) in enumerate(new_params, 1):
        db_storage_key = (p, cyclo_order, tuple(jordan_blocks), D_K)
        kp = kronecker_symbol(D_K, p)
        local_results_map = {}
        max_power_param = len(jordan_blocks) - 1

        if kp == 1: # Split Case
            density_sums_per_block = []
            for exponent, multiplicity in jordan_blocks:
                complete_order = cyclo_order * (p**exponent)
                gl_key = (p, complete_order, multiplicity, max_power_param)
                
                if gl_key not in gl_latt_db:
                    gl_latt_db[gl_key] = enum_with_inv_gl(p, complete_order, multiplicity, max_power_param)
                
                lattices_by_partition = gl_latt_db[gl_key]
                block_densities = []
                for partition in lattices_by_partition:
                    single_block_args = [(p, cyclo_order, ((exponent, multiplicity),), (partition,))]
                    for _, val in orb_int_gl_engine(single_block_args):
                        block_densities.append(val)
                density_sums_per_block.append(block_densities)

            grand_total_density = prod([sum(block) for block in density_sums_per_block])
            if grand_total_density != 0:
                dummy_inv_key = tuple(Mod(0, 2) for _ in jordan_blocks)
                local_results_map[dummy_inv_key] = grand_total_density

        else: # Inert or Ramified Case
            # 代数データのセットアップ
            alg_data = setup_local_algebraic_data(p, cyclo_order, D_K)
            
            partition_options_per_block = []
            for exponent, multiplicity in jordan_blocks:
                complete_order = cyclo_order * (p**exponent)
                hermi_key = (p, complete_order, multiplicity, max_power_param, D_K)
                
                if hermi_key not in hermi_latt_db:
                    hermi_latt_db[hermi_key] = enum_with_inv_hermi(p, complete_order, multiplicity, max_power_param, D_K, alg_data)
                
                partition_options_per_block.append(list(hermi_latt_db[hermi_key].keys()))

            for partition_combination in itertools.product(*partition_options_per_block):
                calculation_args = [(p, cyclo_order, jordan_blocks, partition_combination)]
                for result_args, density_value in orb_int_unitary(calculation_args, hermi_latt_db, D_K):
                    if density_value != 0:
                        inv_key = result_args[3]
                        local_results_map[inv_key] = local_results_map.get(inv_key, 0) + density_value
                        
        orb_int_unitary_db[db_storage_key] = local_results_map

    print(f"Database update complete for D_K = {D_K}.")
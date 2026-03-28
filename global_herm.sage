load("Lfunc.sage")

# ========================================================================
# 1. Helper functions
# ========================================================================
def prodexp(primes_list, exponents):
    res = 1
    for i in range(len(primes_list)):
        res *= primes_list[i]**(exponents[i])
    return res

def totient(primes_list, exponents):
    res = 1
    for i in range(len(primes_list)):
        if exponents[i] > 0:
            res *= (primes_list[i] - 1) * primes_list[i]**(exponents[i] - 1)
    return res

def ppart(m, p):
    k = 0
    x = m
    while x % p == 0:
        x //= p
        k += 1
    return x, k

# ========================================================================
# 2. Generation of conjugate classes and organization of local data
# ========================================================================
def conjclasses_unitary(n, primes_list):
    N = 2 * n
    exponents = [0] * len(primes_list)
    degrees = []
    
    while totient(primes_list, exponents) <= N:
        m = prodexp(primes_list, exponents)
        if m > 0:
            degrees.append([m, totient(primes_list, exponents)])
        exponents[0] += 1
        j = 1
        while j < len(primes_list) and totient(primes_list, exponents) > N:
            exponents[j-1] = 0
            exponents[j] += 1
            j += 1
            
    degrees.sort(key=lambda x: x[0])
    
    def find_partitions(target_dim, idx):
        if target_dim == 0: return [[]]
        if idx >= len(degrees): return []
        
        m, deg = degrees[idx]
        res = []
        for mult in range(target_dim // deg + 1):
            for rest in find_partitions(target_dim - mult * deg, idx + 1):
                if mult > 0:
                    res.append([[m, mult]] + rest)
                else:
                    res.append(rest)
        return res

    return find_partitions(N, 0)

def tot_ram_sub_conj_classes_unitary(cc, p):
    res_dict = {}
    for i, x in enumerate(cc):
        m, k = ppart(x[0], p)
        mult = x[1]
        if m in res_dict:
            res_dict[m].append([k, mult, i])
        else:
            res_dict[m] = [[k, mult, i]]
    res = []
    for m, d in res_dict.items():
        m_red = m // 2 if m % 4 == 2 else m
        res.append([m_red, tuple((x[0], x[1]) for x in d), [x[2] for x in d]])
    return res

def append_poss_local_invs_unitary(poss_local_invs, p, sub_cc, database, D_K):
    p_invs = []
    sub_invs = []

    for s in sub_cc:
        m_red, dims_tup, indices = s
        
        db_key = (p, m_red, tuple(dims_tup), D_K)
        
        if db_key not in database:
            print(f"Warning: Cache miss for {db_key}. This conjugate class drops to Mass 0.")
            return False
            
        orb_dict = database[db_key]
        curr_sub_invs = []
        at_least_one = False
        
        for invs_tup, orb_int in orb_dict.items():
            HW_list = [invs_tup[i] for i in range(len(indices))]
            curr_sub_invs.append([HW_list, orb_int])
            at_least_one = True
            
        if not at_least_one:
            return False
        sub_invs.append(curr_sub_invs)
        
    for L in itertools.product(*sub_invs):
        HW_list_ordered = [None] * sum([len(s[1]) for s in sub_cc])
        orb_int = 1
        for i, s in enumerate(sub_cc):
            sub_HW_list, sub_orb_int = L[i]
            orb_int *= sub_orb_int
            for j, ind in enumerate(s[2]):
                HW_list_ordered[ind] = sub_HW_list[j]
        p_invs.append([HW_list_ordered, orb_int])
        
    poss_local_invs.append(p_invs)
    return True

# ========================================================================
# 3. Calculation of mass
# ========================================================================
def mass_global_term_unitary(cc, D_K):
    res = 1
    for m, n in cc:
        orig_m = m
        if m % 4 == 2: 
            m //= 2
        res *= 2
        if n % 4 == 2:
            res *= (-1)**(euler_phi(m) // 2)
        for l in range(1, n + 1):
            if l % 2 == 0:
                res *= zetarealcyclo_unitary(m, l // 2, D_K)
            else:
                res *= Lfuncquadcyclo_unitary(m, (l - 1) // 2, D_K)
    return res

def mass_unitary(n, primes_list, negdim_goal, cc, database, D_K):
    poss_local_invs = []
    for p in primes_list:
        sub_cc = tot_ram_sub_conj_classes_unitary(cc, p)
        if not append_poss_local_invs_unitary(poss_local_invs, p, sub_cc, database, D_K):
            return 0
    mass = 0
    num_primes = len(primes_list)
    
    for finite_invs in itertools.product(*poss_local_invs):
        orb_int = prod([finite_invs[i][1] for i in range(num_primes)])

        negdim_list = []
        for j, m_and_mult in enumerate(cc):
            m, mult = m_and_mult
            HW_sum = sum([int(finite_invs[i][0][j]) for i in range(num_primes)]) % 2
            poss_negdims = []
            num_real_places = int(euler_phi(m) // 2) if m > 2 else 1
            for archi_negdims in itertools.product(range(mult + 1), repeat=num_real_places):
                if sum(archi_negdims) % 2 == HW_sum:
                    arch_vol = prod([(-1)**(negdim * (mult - negdim)) * binomial(mult, negdim) / (2**mult) 
                                     for negdim in archi_negdims])
                    poss_negdims.append([2 * sum(archi_negdims) if m > 2 else sum(archi_negdims), arch_vol])
            negdim_list.append(poss_negdims)
            
        for archi_negdims in itertools.product(*negdim_list):
            if sum([x[0] for x in archi_negdims]) == negdim_goal:
                mass += orb_int * prod([x[1] for x in archi_negdims])

    global_factor = mass_global_term_unitary(cc, D_K)
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
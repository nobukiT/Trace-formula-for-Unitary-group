load("gllatticeclass.sage") 

def enumsplit(p, m, mult, maxpower):
    """
    GL(N, Q_p) における Z_p-格子（不変部分空間）を全列挙する。
    """
    m2 = m
    k = 0
    result = []
    
    while m2 % p == 0:
        k += 1
        m2 //= p
        
    if k == 0:
        d = 0
        e = 1
    else:
        d = ((p - 1) * k - 1) * p**(k - 1)
        e = p**(k - 1) * (p - 1)

    def generate_I_lists(target_mult, max_idx):
        if max_idx == 0: return [[target_mult]]
        res = []
        for count in range(target_mult + 1):
            for rest in generate_I_lists(target_mult - count, max_idx - 1):
                res.append(rest + [count])
        return res

    ftypes = generate_I_lists(mult, e * maxpower)
    
    for formal_type in ftypes:
        latt = gl_lattice_unr_from_type(p, k, m2, d, e, mult, formal_type)
        result.append(latt)
        
    return result
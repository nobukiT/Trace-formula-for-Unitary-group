load("gllatticeclass.sage") 

def enum_split(p, m, mult, maxpower, ev, qv):
    """
    Enumerate all invariant F_v-lattices in GL(N, F_v) for a given local field F_v.
    
    INPUT:
        p        : a prime number
        m        : complete cyclotomic order
        mult     : dimension or multiplicity
        maxpower : maximum p-power depth parameter
        ev      : ramification index of the local field F_v over Q_p (default: 1)
        qv      : cardinality of the residue field of F_v (default: p)

    OUTPUT:
        A list of `[disc_dummy, lattice_obj]` lists, representing all valid 
        invariant lattices for the given formal types.
    """

    m2 = m
    k = 0
    result = []
    
    while m2 % p == 0:
        k += 1
        m2 //= p
        
    if k == 0:
        d = 0
        absolute_e = 1
    else:
        d = ((p - 1) * k - 1) * p**(k - 1)
        absolute_e = p**(k - 1) * (p - 1)

    relative_e = absolute_e // gcd(absolute_e, ev)

    def generate_lists(target_mult, max_idx):
        if max_idx == 0: return [[target_mult]]
        res = []
        for count in range(target_mult + 1):
            for rest in generate_lists(target_mult - count, max_idx - 1):
                res.append(rest + [count])
        return res

    # The number of slots depends on the relative ramification depth
    num_slots = relative_e * maxpower
    ftypes = generate_lists(mult, num_slots)
    
    for formal_type in ftypes:
        # Pass qv to ensure the local mass term is computed correctly for F_v
        latt = build_gl_lattice(p, k, m2, d, relative_e, mult, formal_type, qv)
        result.append(latt)
        
    return result
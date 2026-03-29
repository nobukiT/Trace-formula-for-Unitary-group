def block_fields_unitary(m, D_K, verbose=False):
    """
    construct the corresponding fields E_i/F_i for the centralizer block.
    OUTPUT:
        dict with keys:
            "m"          : m
            "E_base"     : the quadratic field E = Q(sqrt(D_K))
            "cyclo"      : Q(zeta_m)
            "Ei"         : E(zeta_m)
            "Fi"         : maximal totally real subfield (when Ei is CM)
            "embed_F_to_E": embedding from Fi to Ei
            "split"      : whether Ei is not a field-CM block in the expected sense
            "deg_Ei"     : [Ei:Q]
            "deg_Fi"     : [Fi:Q]
    """
    QQx = QQ['x']
    x = QQx.gen()
    E0 = QuadraticField(D_K, 'w')

    # m=1,2 are degenerate but still valid
    if m == 1:
        L = QQ
    else:
        L = CyclotomicField(m)

    # E(zeta_m)
    Ei, emb1, emb2 = E0.composite_fields(L, both_maps=True)[0]

    out = {
        "m": m,
        "E_base": E0,
        "cyclo": L,
        "Ei": Ei,
        "deg_Ei": Ei.degree(),
    }
    
    Fi, embed_F_to_E = Ei.maximal_totally_real_subfield()
    
    out["Fi"] = Fi
    out["embed_F_to_E"] = embed_F_to_E
    out["deg_Fi"] = Fi.degree()
    out["split"] = False

    return out

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


def is_elliptic(m, D_K):
    """
    Return True iff Phi_m(x), viewed over K[x], has NO GL-type obstruction.
    That is: every irreducible factor f of Phi_m over K is paired with
    inverse_root_polynomial(f) with the same multiplicity.
    Since Phi_m is square-free, this means:
      for every irreducible factor f, inverse_root_polynomial(f) is also
      a factor of Phi_m over K.
    This is the correct criterion for "Phi_m contributes no GL-type factor".
    """
    K = QuadraticField(D_K, 'w')
    R.<x> = PolynomialRing(K)
    phi = R(cyclotomic_polynomial(m))
    fac = factor_dict(phi)

    for f in fac:
        f_inv = inverse_root_polynomial(f)
        if f_inv not in fac:
            return False
    return True

#============================================================================

def prodexp(primes_list, exponents):
    """
    Return prod p_i^e_i
    """
    res = 1
    for i in range(len(primes_list)):
        res *= primes_list[i]**(exponents[i])
    return res

def totient(primes_list, exponents):
    """
    Return phi(prod p_i^e_i)
    """
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

#============================================================================

def conjclasses_unitary(n, primes_list, D_K, verbose=False):
    """
    Enumerate possible Q-characteristic polynomials prod Phi_m(x)^{a_m}
    for finite-order elliptic semisimple stable eigenvalue data
    in U(n,n) over the imaginary quadratic field K.
    OUTPUT FORMAT:
        list of lists
        [
          [[m1,a1],[m2,a2],...],
          [[m1,b1],[m2,b2],...],
          ...
        ]
    Here the total degree is n, i.e.
        sum a_m * phi(m) = n.
    Only those m are allowed for which Phi_m contributes no GL-type factor over K.
    """
    exponents = [0] * len(primes_list)
    degrees = []

    # enumerate m with phi(m) <= n
    while totient(primes_list, exponents) <= n:
        m = prodexp(primes_list, exponents)
        phi_m = totient(primes_list, exponents)

        if m > 0 and phi_m > 0:
            if is_elliptic(m, D_K):
                degrees.append([m, phi_m])
                if verbose:
                    print(f"allowed: m={m}, phi(m)={phi_m}")
            elif verbose:
                print(f"rejected: m={m}, phi(m)={phi_m}")

        exponents[0] += 1
        j = 1
        while j < len(primes_list) and totient(primes_list, exponents) > n:
            exponents[j-1] = 0
            exponents[j] += 1
            j += 1

    degrees.sort(key=lambda x: x[0])

    def find_partitions(target_dim, idx):
        if target_dim == 0:
            return [[]]
        if idx >= len(degrees):
            return []

        m, deg = degrees[idx]
        res = []
        for mult in range(target_dim // deg + 1):
            for rest in find_partitions(target_dim - mult * deg, idx + 1):
                if mult > 0:
                    res.append([[m, mult]] + rest)
                else:
                    res.append(rest)
        return res

    return find_partitions(n, 0)

def tot_ram_sub_conj_classes_unitary(cc, p):
    """
    Input:
        cc = [[m1,a1],[m2,a2],...]
            representing prod Phi_m(x)^a

    Group the cyclotomic factors by their p-prime part m',
    where m = m' * p^k.

    OUTPUT:
        A list of blocks:
        [
          [m_red, ((k1,a1),(k2,a2),...), [original_indices]],
          ...
        ]

    where:
      - m_red is the reduced tame part
      - (ki, ai) records p-adic exponent and multiplicity
      - original_indices remembers where each factor sat in cc
    """
    grouped = {}

    for idx, entry in enumerate(cc):
        m = entry[0]
        mult = entry[1]

        m_prime, k = ppart(m, p)

        if m_prime not in grouped:
            grouped[m_prime] = []

        grouped[m_prime].append([k, mult, idx])

    result = []

    for m_prime, data in grouped.items():
        m_red = m_prime // 2 if m_prime % 4 == 2 else m_prime
        result.append([
            m_red,
            tuple((x[0], x[1]) for x in data),
            [x[2] for x in data]
        ])

    return result

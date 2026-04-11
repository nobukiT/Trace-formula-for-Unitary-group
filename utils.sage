def block_fields_unitary(m, D_K, verbose=False):
    QQx = QQ['x']
    x = QQx.gen()
    E0 = QuadraticField(D_K, 'w')

    if m in (1, 2):
        L = QQ
        Ei = E0
    else:
        L = CyclotomicField(m)
        Ei = E0.composite_fields(L)[0]

    Fi, embed_F_to_E = Ei.maximal_totally_real_subfield()
    
    out = {
        "m": m,
        "E_base": E0,
        "cyclo": L,
        "Ei": Ei,
        "Fi": Fi,
        "embed_F_to_E": embed_F_to_E,
        "deg_Ei": Ei.degree(),
        "deg_Fi": Fi.degree(),
        "split": False
    }

    return out

def inverse_root_polynomial(f):
    R = f.parent()
    coeffs = f.list()
    return R(coeffs[::-1]).monic()

def factor_dict(poly):
    return {g.monic(): e for g, e in poly.factor()}

def is_elliptic(m, D_K):
    K = QuadraticField(D_K, 'w')
    R = PolynomialRing(K, 'x')
    phi = R(cyclotomic_polynomial(m))
    fac = factor_dict(phi)
    for f in fac:
        f_inv = inverse_root_polynomial(f)
        if f_inv not in fac:
            return False
    return True

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

def conjclasses_unitary(n, primes_list, D_K, verbose=False):
    exponents = [0] * len(primes_list)
    degrees = []
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


def exactify_to_QQ(x, max_denominator=10^8, tol=1e-10, prec=200, require_stability=True):
    """
    Convert x to QQ as safely as possible.

    Strategy:
      1. Try exact coercions first.
      2. Try symbolic rational simplification if available.
      3. Only then evaluate numerically at high precision.
      4. Recover a rational only if it is genuinely close and, optionally, stable.
    """
    # 1. Exact coercion first.
    try:
        q = QQ(x)
        if q.denominator() <= max_denominator:
            return q
        raise TypeError(
            "exact rational {} has denominator {} > max_denominator={}".format(
                q, q.denominator(), max_denominator
            )
        )
    except (TypeError, ValueError):
        pass

    # 2. Symbolic simplification, if available.
    if hasattr(x, "simplify_rational"):
        try:
            y = x.simplify_rational()
            q = QQ(y)
            if q.denominator() <= max_denominator:
                return q
        except (TypeError, ValueError):
            pass

    if hasattr(x, "full_simplify"):
        try:
            y = x.full_simplify()
            q = QQ(y)
            if q.denominator() <= max_denominator:
                return q
        except (TypeError, ValueError):
            pass

    # 3. High-precision numeric evaluation.
    RF = RealField(prec)
    try:
        xr = RF(x)
    except (TypeError, ValueError):
        try:
            xr = RF(x.n(prec))
        except Exception as e:
            raise TypeError(
                "unable to numerically evaluate {!r} at prec={}".format(x, prec)
            ) from e

    if xr.is_NaN():
        raise ValueError("NaN encountered while exactifying to QQ")
    if xr.is_infinity():
        raise ValueError("non-finite value encountered while exactifying to QQ")

    # 4. Rational reconstruction with denominator bound.
    q = xr.nearby_rational(max_denominator=max_denominator)
    err = abs(xr - RF(q))

    if err >= RF(tol):
        raise TypeError(
            "unable to safely convert {!r} to QQ: best candidate {} has error {} >= tol={}".format(
                x, q, err, tol
            )
        )

    if require_stability:
        RF2 = RealField(2 * prec)
        try:
            xr2 = RF2(x)
        except (TypeError, ValueError):
            try:
                xr2 = RF2(x.n(2 * prec))
            except Exception as e:
                raise TypeError(
                    "unable to re-evaluate {!r} at prec={}".format(x, 2 * prec)
                ) from e

        if xr2.is_NaN():
            raise ValueError("NaN encountered while exactifying to QQ at higher precision")
        if xr2.is_infinity():
            raise ValueError("non-finite value encountered while exactifying to QQ at higher precision")

        q2 = xr2.nearby_rational(max_denominator=max_denominator)
        err2 = abs(xr2 - RF2(q))

        if q2 != q:
            raise TypeError(
                "rational recovery unstable for {!r}: got {} at prec={}, but {} at prec={}".format(
                    x, q, prec, q2, 2 * prec
                )
            )

        if err2 >= RF2(tol):
            raise TypeError(
                "candidate {} is not stable at higher precision: error {} >= tol={}".format(
                    q, err2, tol
                )
            )

    return QQ(q)
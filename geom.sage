load("global_herm.sage")
load("orb_int_herm.sage")
load("utils.sage")

from itertools import permutations

def generate_Xi(n, a):
    """
    Generate the set Xi_a of minimal-length representatives for W_{M_a} \\ W_G,
    identified with permutations sigma in S_n satisfying

        sigma(a+1) < sigma(a+2) < ... < sigma(n-a).

    Here sigma is represented as a tuple of length n with values in {1,...,n}.
    """
    res = []
    for sigma in permutations(range(1, n + 1)):
        middle = [sigma[i] for i in range(a, n - a)]
        if middle == sorted(middle):
            res.append(sigma)
    return res


def calc_weights(sigma, a, n, k_weights):
    """
    For sigma in Xi_a, compute the weights mu_sigma and nu_{i,sigma}
    appearing in the decomposition

        M_a ~= (Res_{E/Q} GL_1)^a x U(V_a).

    INPUT:
        sigma     : tuple in Xi_a
        a         : integer
        n         : total rank n = p + q
        k_weights : list [k_1, ..., k_n] in descending order

    OUTPUT:
        (mu, nu)
        mu : tuple of length n-2a
        nu : list of pairs (z_{i,1}, z_{i,2}), 1 <= i <= a
    """
    mu = []
    for t in range(1, n - 2 * a + 1):
        pos = sigma[a + t - 1]              # this is sigma(a+t) in 1-based notation
        m_t = k_weights[pos - 1] + a + t - pos
        mu.append(m_t)

    nu = []
    for i in range(1, a + 1):
        pos1 = sigma[i - 1]                 # sigma(i)
        pos2 = sigma[n - i]                 # sigma(n+1-i)

        z_i_1 = k_weights[pos1 - 1] + QQ(n + 1 - 2 * pos1) / QQ(2)
        z_i_2 = k_weights[pos2 - 1] + QQ(n + 1 - 2 * pos2) / QQ(2)
        nu.append((z_i_1, z_i_2))

    return tuple(mu), nu


def cc_eigenvalues(cc, prec=3000):
    """
    Return the eigenvalues over ComplexField(prec) for a stable conjugacy class cc.

    INPUT:
        cc : list such as [[m1, mult1], [m2, mult2], ...]
    """
    CC = ComplexField(prec)
    eigenvals = []
    for m, mult in cc:
        for k in range(1, m + 1):
            if gcd(k, m) == 1:
                root = CC(e ** (2 * I * pi * k / m))
                for _ in range(mult):
                    eigenvals.append(root)
    return eigenvals


def fix_weight_and_get_sign(mu):
    """
    Convert a weight mu to dominant form by Weyl action and return the sign.

    OUTPUT:
        (new_mu, sgn)
        If mu is singular, return (None, 0).
    """
    r = len(mu)
    if r == 0:
        return [], 1

    rho = [QQ(r - 2 * i - 1) / QQ(2) for i in range(r)]
    mu_rho = [mu[i] + rho[i] for i in range(r)]

    if len(set(mu_rho)) < r:
        return None, 0

    sorted_mu_rho = sorted(list(enumerate(mu_rho)), key=lambda x: x[1], reverse=True)

    inv_count = 0
    for i in range(r):
        for j in range(i + 1, r):
            if sorted_mu_rho[i][0] > sorted_mu_rho[j][0]:
                inv_count += 1
    sgn = (-1) ** inv_count

    new_mu = [sorted_mu_rho[i][1] - rho[i] for i in range(r)]
    return new_mu, sgn


def trace_gamma_V_mu(cc, mu, prec=3000):
    """
    Compute tr(gamma | V_mu) using the Weyl character formula.
    """
    if len(mu) == 0:
        return QQ(1)

    fixed_mu, sgn = fix_weight_and_get_sign(mu)
    if sgn == 0:
        return QQ(0)

    eigenvals = cc_eigenvalues(cc, prec)
    if not eigenvals:
        return QQ(sgn)

    min_m = min(fixed_mu)
    shift = 0
    if min_m < 0:
        shift = -min_m

    shifted_mu = [ZZ(x + shift) for x in fixed_mu]
    part = [x for x in shifted_mu if x > 0]

    from sage.combinat.sf.sf import SymmetricFunctions
    CC = ComplexField(prec)
    Sym = SymmetricFunctions(CC)
    s = Sym.schur()

    n_vars = len(eigenvals)
    if n_vars == 0:
        val = CC(1) if not part else CC(0)
    else:
        poly = s(part).expand(n_vars)
        val = poly(eigenvals)

    det_val = prod(eigenvals)
    val *= det_val ** (-shift)

    # In the present setup this should be rational.
    return exactify_to_QQ(sgn * val.real(), prec=max(prec, 300))


def t_ell_unitary(n_a, D_K, primes_list, database, imposed_negdim, mu, prec=3000, verbose=False):
    """
    Compute T_ell(U(V_a), mu) by summing over the mass list for U(V_a).
    """
    if n_a == 0:
        return exactify_to_QQ(trace_gamma_V_mu([], mu, prec=prec), prec=max(prec, 300))

    mass_list = mass_list_unitary(
        n_a,
        primes_list,
        database,
        D_K,
        imposed_negdim=imposed_negdim,
        prec=prec,
        verbose=verbose
    )

    res = QQ(0)
    for cc, mass in mass_list:
        tr = trace_gamma_V_mu(cc, mu, prec=prec)
        res += exactify_to_QQ(mass, prec=max(prec, 300)) * tr

    return res


def t_ell_res_gl1_via_u1(nu_i, D_K, primes_list, database, prec=3000, verbose=False, cache=None):
    """
    Compute T_ell(Res_{E/Q} GL_1, nu_i) via the identification

        E^x \\ A_E^x / A^x  ~=  U_1(Q) \\ U_1(A),

    as stated in the paper.

    If nu_i = (z1, z2), then on U_1 the corresponding algebraic character is
    x |-> x^(z1-z2), so the relevant highest weight for U(1) is mu = (z1-z2).
    """
    z1, z2 = nu_i
    wt = z1 - z2

    try:
        wtZ = ZZ(wt)
    except TypeError:
        # If the weight is not integral, the algebraic character does not descend correctly.
        return QQ(0)

    key = (D_K, tuple(primes_list), wtZ, prec)
    if cache is not None and key in cache:
        return cache[key]

    val = t_ell_unitary(
        n_a=1,
        D_K=D_K,
        primes_list=primes_list,
        database=database,
        imposed_negdim=0,   # U(1) is anisotropic over R
        mu=(wtZ,),
        prec=prec,
        verbose=verbose
    )

    val = exactify_to_QQ(val, prec=max(prec, 300))

    if cache is not None:
        cache[key] = val

    return val


def permutation_sign(sigma):
    """
    Sign epsilon(sigma) of a permutation sigma written as a tuple.
    """
    n = len(sigma)
    inv_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if sigma[i] > sigma[j]:
                inv_count += 1
    return (-1) ** inv_count


def T_geom(p, q, k_weights, D_K, primes_list, database, prec=3000, verbose=False):
    """
    Compute the geometric side contribution according to Proposition
    'parabolic_contribution_Ma' in the paper:

        sum_a
          [ (-1)^{a + (q-a)(q-a+1)/2} 2^{q-2a} (q-a)! / a! ]
          sum_{sigma in Xi_a}
            epsilon(sigma)
            T_ell(U(V_a), mu_sigma)
            prod_i T_ell(Res_{E/Q} GL_1, nu_{i,sigma}).

    INPUT:
        p, q       : signature parameters, n = p + q
        k_weights  : highest weight (k_1, ..., k_n), in descending order
        D_K        : discriminant of the imaginary quadratic field E
        primes_list: list of primes used in the local computation
        database   : orbital integral database
        prec       : precision for character computations
        verbose    : boolean
    """
    n = p + q
    total_T = QQ(0)
    gl1_cache = {}

    for a in range(min(p, q) + 1):
        n_a = n - 2 * a
        q_a = q - a
        Xi_a = generate_Xi(n, a)

        sum_sigma = QQ(0)

        for sigma in Xi_a:
            eps_sigma = permutation_sign(sigma)
            mu_sigma, nu_sigma = calc_weights(sigma, a, n, k_weights)

            gl1_prod = QQ(1)
            for nu_i in nu_sigma:
                val_i = t_ell_res_gl1_via_u1(
                    nu_i,
                    D_K=D_K,
                    primes_list=primes_list,
                    database=database,
                    prec=prec,
                    verbose=verbose,
                    cache=gl1_cache
                )
                if val_i == 0:
                    gl1_prod = QQ(0)
                    break
                gl1_prod *= val_i

            if gl1_prod == 0:
                continue

            t_unitary = t_ell_unitary(
                n_a=n_a,
                D_K=D_K,
                primes_list=primes_list,
                database=database,
                imposed_negdim=q_a,
                mu=mu_sigma,
                prec=prec,
                verbose=verbose
            )

            sum_sigma += eps_sigma * t_unitary * gl1_prod

        coef = ((-1) ** (a + (q - a) * (q - a + 1) // 2)) * (2 ** (q - 2 * a))  * factorial(q - a) / factorial(a)

        total_T += coef * sum_sigma

        if verbose:
            print("a =", a, "coef =", coef, "partial contribution =", exactify_to_QQ(coef * sum_sigma))

    return exactify_to_QQ(total_T, prec=max(prec, 300))
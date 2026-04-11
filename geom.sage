load("global_herm.sage")
load("orb_int_herm.sage")
load("utils.sage")
from itertools import permutations
def get_E_class_number_and_roots(D_K):
    """
    Computes the class number h_E and the number of roots of unity w_E 
    for the imaginary quadratic field E = Q(sqrt(D_K)).

    INPUT:
        D_K : integer; the discriminant of the imaginary quadratic field.
    
    OUTPUT:
        tuple (h_E, w_E); 
        - h_E : integer; class number of E.
        - w_E : integer; number of roots of unity in E.
    """
    E.<a> = QuadraticField(D_K)
    h_E = E.class_number()
    w_E = len(E.roots_of_unity())
    return h_E, w_E

def generate_Xi(n, a):
    """
    Generates the set of minimal-length representatives Xi_a for the quotient 
    W_{M_a} \ W_G. These are permutations sigma in S_n satisfying 
    sigma(a+1) < sigma(a+2) < ... < sigma(n-a).

    INPUT:
        n : integer; the total dimension n = p + q.
        a : integer; the number of GL_1 factors in the Levi subgroup M_a.

    OUTPUT:
        list of permutations (tuples) representing the elements of Xi_a.
    """
    res = []
    for sigma in permutations(range(1, n+1)):
        # Check if the elements from index 'a' to 'n-a-1' are in ascending order (0-based)
        sub = [sigma[i] for i in range(a, n-a)]
        if list(sub) == sorted(sub):
            res.append(sigma)
    return res

def calc_weights(sigma, a, n, k_weights):
    """
    Calculates the highest weights mu_sigma (for the U(V_a) component) 
    and nu_i_sigma (for the GL_1 components) based on the permutation sigma 
    and the overall highest weight k_weights.

    INPUT:
        sigma     : tuple; a permutation from Xi_a.
        a         : integer; the number of GL_1 factors.
        n         : integer; the total dimension.
        k_weights : list; the highest weight lambda = (k_1, ..., k_n).

    OUTPUT:
        tuple (mu, nu);
        - mu : tuple; the highest weight mu_sigma for the U(V_a) factor.
        - nu : list of tuples; the highest weights nu_{i, sigma} for the GL_1 factors.
    """
    # mu_sigma (Weight for the U(V_a) component)
    mu = []
    for t in range(1, n - 2*a + 1):
        idx = sigma[a + t - 1] - 1  # 0-based index
        m_t = k_weights[idx] + a + t - sigma[a + t - 1]
        mu.append(m_t)
    
    # nu_sigma (Weights for the GL_1 components)
    nu = []
    for i in range(1, a + 1):
        idx1 = sigma[i - 1] - 1
        idx2 = sigma[n + 1 - i - 1] - 1
        z_i_1 = k_weights[idx1] + (n + 1 - 2*sigma[i - 1]) / 2
        z_i_2 = k_weights[idx2] + (n + 1 - 2*sigma[n + 1 - i - 1]) / 2
        nu.append((z_i_1, z_i_2))
        
    return tuple(mu), nu

def t_ell_gl1(nu_i, h_E, w_E, num_ram):
    """
    Computes the elliptic contribution T_{ell} for the Res_{E/Q} GL_1 factor.
    INPUT:
        nu_i    : tuple (z1, z2); the highest weight for the GL_1 factor.
        h_E     : integer; class number of E.
        w_E     : integer; number of roots of unity in E.
        num_ram : integer; number of ramified primes in E/Q.

    OUTPUT:
        rational;
        the elliptic contribution of the GL_1 factor, adjusted for orbital integrals.
    """
    z1, z2 = nu_i
    if (z1 - z2) % w_E == 0:
        return QQ(2 * h_E) / QQ(2**num_ram)
    else:
        return QQ(0)

def cc_eigenvalues(cc, prec=200):
    """
    Returns the eigenvalues over the complex field for a given conjugacy class cc.

    INPUT:
        cc   : list; the stable conjugacy class represented as [[m1, mult1], [m2, mult2], ...].
        prec : integer; precision for the complex numbers.

    OUTPUT:
        list of complex numbers; the eigenvalues corresponding to the conjugacy class.
    """
    CC = ComplexField(prec)
    eigenvals = []
    for m, mult in cc:
        for k in range(1, m+1):
            if gcd(k, m) == 1:
                root = CC(e^(2 * I * pi * k / m))
                for _ in range(mult):
                    eigenvals.append(root)
    return eigenvals

def fix_weight_and_get_sign(mu):
    """
    Converts a weight mu to a dominant weight (descending order) and returns 
    the corresponding sign character from the Weyl group action.

    INPUT:
        mu : tuple or list; the weight to be processed.

    OUTPUT:
        tuple (new_mu, sgn);
        - new_mu : list; the dominant weight.
        - sgn    : integer; the sign character (1, -1, or 0 if singular).
    """
    r = len(mu)
    if r == 0:
        return [], 1
    
    rho = [(r - 2*i - 1)/2 for i in range(r)]
    mu_rho = [mu[i] + rho[i] for i in range(r)]
    
    if len(set(mu_rho)) < r:
        return None, 0
        
    sorted_mu_rho = sorted(list(enumerate(mu_rho)), key=lambda x: x[1], reverse=True)
    
    inv_count = 0
    for i in range(r):
        for j in range(i+1, r):
            if sorted_mu_rho[i][0] > sorted_mu_rho[j][0]:
                inv_count += 1
    sgn = (-1)**inv_count
    
    new_mu = [sorted_mu_rho[i][1] - rho[i] for i in range(r)]
    return new_mu, sgn

def trace_gamma_V_mu(cc, mu, prec=200):
    """
    Computes the trace of a finite-dimensional representation V_mu 
    evaluated at a conjugacy class gamma using the Weyl character formula.

    INPUT:
        cc   : list; the stable conjugacy class [[m, mult], ...].
        mu   : tuple or list; the highest weight of the representation.
        prec : integer; precision for the complex field.

    OUTPUT:
        exact rational number; the value of the trace evaluated at gamma.
    """
    if len(mu) == 0:
        return QQ(1)
        
    fixed_mu, sgn = fix_weight_and_get_sign(mu)
    if sgn == 0:
        return QQ(0)
        
    eigenvals = cc_eigenvalues(cc, prec)
    if not eigenvals:
        return exactify_to_QQ(sgn)
        
    min_m = min(fixed_mu)
    shift = 0
    if min_m < 0:
        shift = -min_m
        
    shifted_mu = [int(x + shift) for x in fixed_mu]
    part = [x for x in shifted_mu if x > 0]  # Filter out trailing zeros for Partition
    
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
    return exactify_to_QQ(sgn * val.real())
    
def t_ell_U_Va(n_a, D_K, primes_list, database, imposed_negdim, mu, prec=200, verbose=False):
    """
    Computes the elliptic contribution T_{ell} for the smaller unitary group U(V_a) 
    using the pre-existing mass_list_unitary function.

    INPUT:
        n_a            : integer; dimension of the smaller unitary group U(V_a).
        D_K            : integer; discriminant of the imaginary quadratic field.
        primes_list    : list; the primes considered for local data.
        database       : dict; the local orbital integral database.
        imposed_negdim : integer; the signature q_a.
        mu             : tuple; the highest weight for the U(V_a) component.
        prec           : integer; computation precision.
        verbose        : boolean; whether to print progress.

    OUTPUT:
        exact rational number; the elliptic contribution T_{ell}(U(V_a), mu_sigma).
    """
    if n_a == 0:
        return exactify_to_QQ(trace_gamma_V_mu([], mu, prec))
        
    mass_list = mass_list_unitary(n_a, primes_list, database, D_K, imposed_negdim=imposed_negdim, prec=prec, verbose=verbose)
    
    res = QQ(0)
    for cc, mass in mass_list:
        tr = trace_gamma_V_mu(cc, mu, prec=prec)
        res += mass * tr
    return res

def T_geom(p, q, k_weights, D_K, primes_list, database, prec=200, verbose=False):
    """
    Computes the geometric side T_{geom} of the Arthur-Selberg trace formula 
    for the unitary group U(p,q) as described in the paper.

    INPUT:
        p           : integer; positive signature of the Hermitian form.
        q           : integer; negative signature of the Hermitian form.
        k_weights   : list; highest weight lambda = (k_1, ..., k_n) in descending order.
        D_K         : integer; discriminant of the imaginary quadratic field E.
        primes_list : list; the primes to be checked.
        database    : dict; the precomputed local orbital integral database.
        prec        : integer; precision for complex numbers.
        verbose     : boolean; whether to output debug information.

    OUTPUT:
        exact rational number; the total evaluated value of the geometric side T_{geom}.
    """
    n = p + q
        
    h_E, w_E = get_E_class_number_and_roots(D_K)
    num_ram = len(prime_divisors(abs(D_K)))
    total_T = QQ(0)
    
    for a in range(min(p, q) + 1):
        n_a = n - 2*a
        q_a = q - a
        
        sum_sigma = QQ(0)
        Xi_a = generate_Xi(n, a)
        
        for sigma in Xi_a:
            # Compute the sign character epsilon(sigma) of the permutation
            inv_count = 0
            for i in range(n):
                for j in range(i+1, n):
                    if sigma[i] > sigma[j]:
                        inv_count += 1
            sgn = (-1)**inv_count
            
            mu, nu = calc_weights(sigma, a, n, k_weights)
            
            # Product of GL_1 contributions (skip U(V_a) heavy computation if 0)
            gl1_prod = QQ(1)
            for nu_i in nu:
                val = t_ell_gl1(nu_i, h_E, w_E, num_ram)
                if val == 0:
                    gl1_prod = QQ(0)
                    break
                gl1_prod *= val
            
            if gl1_prod == 0:
                continue
                
            # Contribution from U(V_a)
            t_ell_U = t_ell_U_Va(n_a, D_K, primes_list, database, q_a, mu, prec=prec, verbose=verbose)
            
            sum_sigma += sgn * t_ell_U * gl1_prod
            
        coef = ((-1)**(a + (q-a)*(q-a+1)//2)) * QQ(2**(q-2*a)) * factorial(q-a) / factorial(a)
        total_T += coef * sum_sigma
        
    return exactify_to_QQ(total_T)
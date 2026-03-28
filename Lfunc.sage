# m>=3, m != 2 mod 4
# F = QQ(eta_m), k>=1
# computes zeta_F(1-2k)
# note: m=3 and 4 yields the Riemann zeta function
def zetarealcyclo(m, k):
    """
    Q(zeta_m) の最大実部分体のゼータ関数の負の偶数 (-2k+1) での特殊値に対応する因子。
    偶数指標 (chi(-1) == 1) を持つディリクレ指標に対するベルヌーイ数の積を計算する。
    """
    G = DirichletGroup(m)
    res = 1
    for e in G:
        if e(-1) == 1:
            chi_prim = e.primitive_character()
            # -B_{2k, chi} / 2k
            res *= -chi_prim.bernoulli(2 * k) / (2 * k)
    return QQ(res)

def Lfuncquadcyclo(m, k):
    """
    Q(zeta_m) の相対 L-関数の負の奇数 (-2k) での特殊値に対応する因子。
    奇数指標 (chi(-1) == -1) を持つディリクレ指標に対するベルヌーイ数の積を計算する。
    """
    G = DirichletGroup(m)
    res = 1
    for e in G:
        if e(-1) == -1:
            chi_prim = e.primitive_character()
            # -B_{2k+1, chi} / (2k+1)
            res *= -chi_prim.bernoulli(2 * k + 1) / (2 * k + 1)
    return QQ(res)

def zetarealcyclo_unitary(m, k, D_K):
    """
    U(n,n) 用: 偶数次 (l が偶数) における L関数の特殊値因子。
    chi_K^even = 1 なので、通常の zetarealcyclo と同じ（円分体の偶指標の積）。
    """
    G = DirichletGroup(m)
    res = 1
    for e in G:
        if e(-1) == 1:
            chi_prim = e.primitive_character()
            res *= -chi_prim.bernoulli(2 * k) / (2 * k)
    return QQ(res)

def kronecker_dirichlet_character(D, M=None):
    """
    判別式 D に対応する Kronecker 指標 χ_D を
    modulus M 上の Dirichlet character として返す。

    比較は gcd(n,M)=1 の n のみで行う。
    """
    N = abs(D)
    if M is None:
        M = N
    if M % N != 0:
        raise ValueError("Need abs(D) | M")
    
    G = DirichletGroup(M, QQbar)
    
    for chi in G:
        ok = True
        for n in range(M):
            if gcd(n, M) == 1:
                if chi(n) != kronecker(D, n):
                    ok = False
                    break
        if ok:
            return chi
    
    raise ValueError(f"Kronecker character for D={D} not found modulo {M}")


def lift_character_to_modulus(chi, M):
    """
    modulus m の Dirichlet character chi を
    modulus M (m | M) に持ち上げる。
    """
    m = chi.modulus()
    if M % m != 0:
        raise ValueError("Need modulus(chi) | M")
    
    G = DirichletGroup(M, QQbar)
    
    for psi in G:
        ok = True
        for a in range(M):
            if gcd(a, M) == 1:
                if psi(a) != chi(a % m):
                    ok = False
                    break
        if ok:
            return psi
    
    raise ValueError("Could not lift character")


def Lfuncquadcyclo_unitary(m, k, D_K):
    """
    U(n,n) 用: 奇数次の L関数特殊値因子
    """
    M = lcm(m, abs(D_K))
    
    G_m = DirichletGroup(m, QQbar)
    chi_K = kronecker_dirichlet_character(D_K, M)
    
    res = 1
    
    for e in G_m:
        e_lift = lift_character_to_modulus(e, M)
        chi_twist = e_lift * chi_K
        
        if chi_twist(-1) == -1:
            chi_prim = chi_twist.primitive_character()
            res *= -chi_prim.bernoulli(2 * k + 1) / (2 * k + 1)
    
    return QQ(res)

def card_unitary_group(p, m, dim):
    f = 1
    x = Mod(p, m)
    while x != Mod(-1, m):
        f += 1
        x *= p
        if x == Mod(1, m):
            raise ValueError("Extension should be unramified; it is split!")
    if f != euler_phi(m) // 2:
        raise ValueError("Error in card_unitary_group: f is wrong.")
    q = p**f
    return q**(dim**2) * prod([(1 - (-q)**(-l)) for l in range(1, dim+1)])


def card_general_linear_group(p, m, dim):
    f = 1
    x = Mod(p, m)
    while x != Mod(1, m):
        f += 1
        x *= p
        if x == Mod(-1, m):
            raise ValueError("Extension should be split; it is unramified!")
    if f != euler_phi(m) // 2:
        raise ValueError("Error in card_general_linear_group: f is wrong.")
    q = p**f
    return q**(dim**2) * prod([1 - q**(-l) for l in range(1, dim+1)])
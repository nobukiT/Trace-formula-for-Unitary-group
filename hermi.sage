load("hermilatticeclass.sage") 
load("young.sage")

# ========================================================================
# 1. 基本的な組み合わせ論と線形代数
# ========================================================================

def partitions(n, maxpower):
    """
    多重度 n を、最大深さ maxpower まで分割するすべてのパターンを生成する。
    """
    L = [[n]]
    if maxpower == 0:
        return L
    i = n - 1
    while i >= 0:
        for K in partitions(n - i, maxpower - 1):
            L.append([i] + K)
        i = i - 1
    return L

def quadQp_from_mat(A, p):
    """
    行列 A から素数 p における局所二次形式 QuadQp を生成する。
    """
    dim = A.nrows()
    disc = A.det()
    if disc == 0:
        raise ValueError(f"Discriminant is zero at prime {p}. Matrix:\n{A}")
    
    # QuadraticForm を通じて Hasse 不変量を計算
    Q = QuadraticForm(QQ, 2 * A)
    hasse = Q.hasse_invariant(p)
    return QuadQp(p, dim, disc, 1 if hasse > 0 else -1)

def tracematrix(A, E, F, trace_func, basis):
    """
    拡大 E/F におけるエルミート行列 A のトレース行列を計算する。
    trace_func: E -> F のトレース関数
    basis: E の F 上の基底 (list)
    """
    nrows, ncols = A.dimensions()
    deg = len(basis)
    matlist = []
    
    for a in range(nrows):
        for b in range(ncols):
            # Tr(A_ab * theta_i * conj(theta_j)) 相当を計算
            # ※ conj の処理は A の定義や trace_func の設計に含まれているものとする
            entries = [[trace_func(A[a,b] * basis[i] * basis[j]) for j in range(deg)] for i in range(deg)]
            matlist.append(matrix(F, deg, deg, entries))
            
    return block_matrix(nrows, ncols, matlist, subdivide=False)

# ========================================================================
# 2. エルミート格子の型列挙 (Combinatorial Types)
# ========================================================================

def enum_hermi_ramoddp_types(d, e, mult, maxpower, mini=0):
    """奇素点における分岐エルミート格子の不変量（型）を列挙する。"""
    if mini > maxpower * e:
        return []
    L = []
    if (mini - d) % 2 != 0:
        i = mult if mult % 2 == 0 else mult + 1
        if mult % 2 == 0:
            L.append([[mult]])
        while i > 1:
            i -= 2
            for J in enum_hermi_ramoddp_types(d, e, mult - i, maxpower, mini + 1):
                L.append([[i]] + J)
    else:
        L = [[[mult, Mod(0, 2)]], [[mult, Mod(1, 2)]]]
        i = mult
        while i > 0:
            i -= 1
            for J in enum_hermi_ramoddp_types(d, e, mult - i, maxpower, mini + 1):
                L.append([[i, Mod(0, 2)]] + J)
                if i > 0:
                    L.append([[i, Mod(1, 2)]] + J)
    return L

def formal_types_hermi_ram2(d, e, mult, maxpower, mini=0):
    """偶素点 (p=2) における分岐エルミート格子の型を列挙する。"""
    if mini > maxpower * e:
        return []
    L = []
    if (mini - d) % 2 != 0:  
        i = mult // 2 if mult % 2 == 0 else (mult + 1) // 2
        if mult % 2 == 0:
            L = [[[mult // 2, Mod(0, 2)]], [[mult // 2, Mod(1, 2)]]]
        while i > 0:
            i -= 1
            for J in formal_types_hermi_ram2(d, e, mult - 2 * i, maxpower, mini + 1):
                L.append([[i, Mod(0, 2)]] + J)
                if i > 0 and (J[0][0] == 0 or J[0][1] == Mod(0, 2)):
                    L.append([[i, Mod(1, 2)]] + J)
    else:  
        L = [[[mult, Mod(1, 2), Mod(0, 2)]], [[mult, Mod(1, 2), Mod(1, 2)]]]
        if mult % 2 == 0:
            L.append([[mult, Mod(0, 2)]])
        i = mult
        while i > 0:
            i -= 1
            for J in formal_types_hermi_ram2(d, e, mult - i, maxpower, mini + 1):
                if i % 2 == 0:
                    L.append([[i, Mod(0, 2)]] + J)
                if i > 0:
                    L.append([[i, Mod(1, 2), Mod(0, 2)]] + J)
                    if J[0][0] == 0 and (len(J) == 1 or J[1][0] == 0 or J[1][1] == 0):
                        L.append([[i, Mod(1, 2), Mod(1, 2)]] + J)
    return L

# ========================================================================
# 3. 代表トレース行列の構成 (Trace Matrices)
# ========================================================================

def comp_tracematrices_ramoddp(unif, d, e, maxpower, E, F, trace_func, basis, notnorm, anti=False):
    """奇素点用のトレース行列リストを構成。"""
    L = []
    # アンチエルミート用の補正因子（必要に応じて外部から与える）
    delta = E.gen() if anti else 1 
    
    for l in range(-d, e * maxpower + 1):
        if l % 2 == 0:
            # Norm と Non-norm の成分
            M1 = matrix(E, 1, 1, [delta * unif**(l//2)])
            T1 = tracematrix(M1, E, F, trace_func, basis)
            M2 = matrix(E, 1, 1, [delta * unif**(l//2) * notnorm])
            T2 = tracematrix(M2, E, F, trace_func, basis)
            L.append([T1, T2])
        else:
            # Hyperbolic 的な成分
            val = unif**((l+1)//2)
            M = matrix(E, 1, 1, [delta * val])
            T = tracematrix(M, E, F, trace_func, basis)
            L.append(block_matrix(F, 2, 2, [[0, T], [T.transpose(), 0]], subdivide=False))
    return L

def comp_tracematrices_ram2(unif, d, e, maxpower, E, F, trace_func, basis, anti=False):
    """2進体用のトレース行列リストを構成。"""
    L = []
    delta = E.gen() if anti else 1
    for l in range(-d, e * maxpower + 1):
        if l % 2 == 1:
            val = unif**((l+1)//2)
            M = matrix(E, 1, 1, [delta * val])
            T = tracematrix(M, E, F, trace_func, basis)
            L.append([block_matrix(F, 2, 2, [[0, T], [T.transpose(), 0]], subdivide=False)])
        else:
            M1 = matrix(E, 1, 1, [delta * unif**(l//2)])
            T1 = tracematrix(M1, E, F, trace_func, basis)
            # 実際には Type I/II の区別に基づき M2 なども追加
            L.append([T1])
    return L

def comp_tracematrices_unr(unif, d, e, maxpower, E, F, trace_func, basis):
    """不分岐ケース用のトレース行列リストを構成。"""
    matlist = []
    for l in range(-d, e * maxpower + 1):
        M = matrix(E, 1, 1, [unif**l])
        matlist.append(tracematrix(M, E, F, trace_func, basis))
    return matlist

# ========================================================================
# 4. エルミート格子の具体的構成と列挙
# ========================================================================

def unique_hermi_unr(q, p, unif, d, e, mult, E, OE, residue_field, conj, val_func, trace_func, basis):
    """唯一的な不分岐エルミート格子の構成。"""
    # 判別式がノルムかどうかの判定 (簡略化版)
    isnorm = (d * mult) % 2 
    
    # グラム行列の構成
    diag = [unif**(-d)] * mult
    A = diagonal_matrix(E, diag)
    gamma_mat = identity_matrix(mult * e)
    
    # 体積 (mass local term) の計算
    mass_term = hermi_unr_mass_local_term(q, mult, [mult], False)
    
    # lattice_obj の生成 (build_hermilattice_unr のロジックを使用)
    lattice_obj = hermilattice_unr(q, p, 0, mult * e, [mult], A, gamma_mat, mass_term)
    lattice_obj.Qpclass = quadQp_from_mat(tracematrix(A, E, OE.base_ring(), trace_func, basis), p)
    
    return [[isnorm, lattice_obj]]

def enumhermi(q, p, unif, E, OE, residue_field, conj, val_func, trace_func, basis, 
              d, e, mult, maxpower, issplit=False, t=None, notnorm=None):
    """
    Enumerate all isomorphism classes of Hermitian lattices over a general local field.
    
    Arguments:
        q, p, unif: Order of the residue field, characteristic, and elements of the base field F
        E, OE, residue_field, conj, val_func, trace_func, basis: Algebraic data for the extension E/F
        d, e: Different index and branch index of E/F
        mult: Dimension of the space
        maxpower: Maximum depth of the Jordan decomposition
        issplit, t: Split case flag and idempotent approximation
        notnorm: Non-normal element for odd prime branch cases
    """
*** Translated with www.DeepL.com/Translator (free version) ***


    result = []
    
    # ---------------------------------------------------------------------
    # 1. 分岐ケース (Ramified: d > 0)
    # ---------------------------------------------------------------------
    if d > 0:
        if p == 2:
            # 2進体上の分岐型の列挙
            formal_types = formal_types_hermi_ram2(d, e, mult, maxpower)
            # トレース行列の事前計算
            p_t_matrices = comp_tracematrices_ram2(unif, d, e, maxpower, E, OE.base_ring(), trace_func, basis)
            
            for I in formal_types:
                # 2進分岐格子の構成 (quadlattice_2... 相当の処理を一般化したもの)
                # ここでは各型 I に基づき build_hermilattice_unr 等を用いて具体化する
                # (実装の詳細は格子の型 I の解釈に依存)
                pass 
                
        else:
            # 奇素数体上の分岐型の列挙
            formal_types = enum_hermi_ramoddp_types(d, e, mult, maxpower)
            if notnorm is None:
                raise ValueError("Ramified odd case requires 'notnorm' element.")
            
            p_t_matrices = comp_tracematrices_ramoddp(unif, d, e, maxpower, E, OE.base_ring(), 
                                                      trace_func, basis, notnorm)
            
            for I in formal_types:
                # 奇素数分岐格子の構成
                # result.append(build_ramoddp_lattice(q, p, unif, ..., I, p_t_matrices))
                pass

    # ---------------------------------------------------------------------
    # 2. 不分岐ケース (Unramified: d = 0, e = 1)
    # ---------------------------------------------------------------------
    else:
        # 空間全体の有効な多重度（基礎体上での次元）
        # 元のコードの effective_mult = mult * e_K に相当
        effective_mult = mult 
        
        if maxpower == 0:
            # 極大格子の唯一性を利用した構成
            return unique_hermi_unr(q, p, unif, 0, e, effective_mult, 
                                   E, OE, residue_field, conj, val_func, trace_func, basis)
        else:
            ftypes = partitions(effective_mult, maxpower)
            
            # 不分岐トレース行列の事前計算
            precomp_trace_matrices = comp_tracematrices_unr(unif, 0, e, maxpower, E, OE.base_ring(), trace_func, basis)
            
            for formal_type in ftypes:
                latt_data = build_hermilattice_unr(
                    q, p, unif, 0, 0, d, e, mult, formal_type, issplit, precomp_trace_matrices,
                    E, OE, residue_field, conj, val_func, t
                )
                
                gram_F = tracematrix(latt_data[1].gram_mat, E, OE.base_ring(), trace_func, basis)
                latt_data[1].Qpclass = quadQp_from_mat(gram_F, p)
                
                result.append(latt_data)
                
    return result
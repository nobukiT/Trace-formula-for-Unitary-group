load("hermi_density.sage")
load("young.sage")

# q: 基礎体の剰余体の位数 (一般の局所体 F_p の剰余体サイズ)
# p: 標数 (characteristic)
# unif: 基礎体の素元 (ユニフォーマイザ \varpi)
# E, OE, residue_field, conj, val_func, t: ジョルダン分解に渡すための代数的構造
# k>=0
# d is the exponent of the different of F_p(zeta_p^k)/F_p
def build_hermilattice_unr(q, p, unif, k, m2, d, e, dim, I, issplit, precomp_trace_matrices, 
                           E, OE, residue_field, conj, val_func, t=None):
    """
    一般の局所代数体上の不分岐（または完全分岐）エルミート格子のデータと局所体積を計算する。
    """
    diag_gram_mat = []
    
    for i, count in enumerate(I):
        if count > 0:
            # 修正: p**(i-d) ではなく、与えられた素元 unif**(i-d) を使用する
            val = precomp_trace_matrices[i] if k > 0 else unif**(i - d)
            diag_gram_mat.extend([val] * count)
            
    if k > 0:
        gram_mat0 = block_diagonal_matrix(diag_gram_mat)
        
        # 一般化されたジョルダン分解を呼び出す
        newbasis, gram_mat, type_list = jordanhermi_generalized(
            gram_mat0, issplit, E, OE, residue_field, unif, conj, val_func, t
        )
        
        # 注意: companion_matrix は特定の多項式に依存します。
        # 完全に一般の局所体上で計算する場合、この cyclotomic_polynomial は
        # 基礎体上の対応する極小多項式に置き換える必要があります。
        comp_mat = companion_matrix(cyclotomic_polynomial(p**k), format='right')
        gamma_mat0 = block_diagonal_matrix([comp_mat] * dim, subdivide=False)
        gamma_mat = newbasis**(-1) * gamma_mat0 * newbasis
    else:
        # ZZ への強制キャストを削除し、E または OE の元として自然に扱わせる
        gram_mat = diagonal_matrix(diag_gram_mat)
        gamma_mat = identity_matrix(dim)
        type_list = I

    # 1. 局所密度の分母（体積）を計算（一般化された q を渡す）
    mass_local_term = hermi_unr_mass_local_term(q, dim, I, issplit)
    
    # 2. Mod() オブジェクトのループ内生成を排除し、純粋な整数演算で高速化
    if issplit:
        disc_isnorm = 0
    else:
        disc_isnorm = sum((i - d) * count for i, count in enumerate(I)) % 2
        
    # クラスのインスタンス化 (q を渡す)
    lattice_obj = hermilattice_unr(q, p, m2, dim * e, type_list, gram_mat, gamma_mat, mass_local_term)
    
    return [disc_isnorm, lattice_obj]


class hermilattice_unr:
    """
    一般の局所体上の不分岐エルミート格子の情報を保持するデータクラス
    """
    def __init__(self, q, p, m2, dim, type_list, gram_mat, gamma_mat, mass_local_term, compute_reduction=True):
        self.q = q  # 剰余体の位数
        self.p = p  # 標数
        self.m2 = m2
        self.dim = dim
        self.type_list = type_list
        self.gram_mat = gram_mat
        self.gamma_mat = gamma_mat
        self.mass_local_term = mass_local_term
        
        self.is_simple = False

        # 3. 剰余体上でのヤング図形(Jordan標準形)の計算
        if compute_reduction and len(type_list) <= 2:
            self.is_simple = True
            
            if len(type_list) <= 1:
                self.reduction_dim = 0
                self.reduction_young = []
            else:
                self.reduction_dim = type_list[1]
                
                # 【重要な変更点】
                # 有理数体に依存した CyclotomicField(m2) からの素イデアル分解を廃止し、
                # 直接、位数 q の有限体 (剰余体) を生成します。
                # これにより、任意の局所体上で動作するようになります。
                residue_field_local = GF(q, 'a')
                
                # 自己同型行列から部分行列を抽出し、剰余体へ落として -I する
                # （gamma_mat の成分が整数の場合、自動的に有限体へ射影・還元されます）
                sub_gamma = gamma_mat.submatrix(type_list[0], type_list[0])
                reduced_matrix = matrix(residue_field_local, self.reduction_dim, self.reduction_dim, sub_gamma) - 1
                
                self.reduction_young = young_list(reduced_matrix)

    def __str__(self):
        return (f"Hermitian lattice over local field with residue field F_{self.q}, "
                f"rank {self.dim}, type {self.type_list}, "
                f"local term in mass formula: {self.mass_local_term}")

    def __repr__(self):
        return str(self)
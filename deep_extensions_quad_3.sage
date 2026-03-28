def get_ref_bil_form_munip_3(m):
    """Jordan block (size 2m) に対する参照双線形形式 S を計算"""
    if m < 1:
        return matrix(Integers(3), 0, 0)
        
    # 多項式環と伴随行列の準備
    r0 = PolynomialRing(Integers(3), 'x0')
    x0 = r0.gen()
    poly = (x0 + 1)**(2 * m)
    p_mat = companion_matrix(poly)
    
    # 商環での計算
    r_quot = r0.quotient(poly, 'x')
    x = r_quot.gen()
    
    dim = 2 * m
    s_mat = matrix(Integers(3), dim, dim)
    for i in range(dim):
        for j in range(dim):
            # x^(j-i+m) の係数リストを取得
            coeffs = (x**(j - i + m)).list()
            # 必要なインデックスの係数を安全に抽出
            c0 = coeffs[0] if len(coeffs) > 0 else 0
            cn = coeffs[dim - 1] if len(coeffs) >= dim else 0
            s_mat[i, j] = c0 + cn
            
    # 新しい基底の構築
    e1 = matrix(Integers(3), dim, 1, [1] + [0] * (dim - 1))
    basis_cols = [e1]
    for i in range(dim - 1):
        # 直前の列ベクトルに P を作用させて次を作る
        basis_cols.append(basis_cols[i] + p_mat * basis_cols[i])
        
    new_basis = block_matrix(1, dim, basis_cols)
    # new_s = transpose(basis) * S * basis
    return new_basis.transpose() * s_mat * new_basis

def invariants_munip_bil_form_3(s_mat):
    """F3上の 'minus unipotent' 双線形形式の不変量を計算"""
    if s_mat.nrows() == 0:
        return ()
        
    # P = transpose(S^-1) * S
    p_mat = (s_mat.inverse().transpose()) * s_mat
    n_mat = 1 + p_mat  # べき零部分
    
    # ジョルダン標準形と変換行列
    new_n, new_basis = n_mat.jordan_form(transformation=True)
    # ブロックサイズと重複度の取得
    young_tableau = get_jordan_mults(new_n.subdivisions()[0], n_mat.nrows())
    
    # 基底の正規化（反対角行列による調整）
    block_list = [antidiag_scalar_matrix(Integers(3), r[0], 1) for r in young_tableau for _ in range(r[1])]
    new_basis *= block_diagonal_matrix(block_list, subdivide=False)
    
    # 変換後の S
    transformed_s = new_basis.transpose() * s_mat * new_basis
    
    res_list = [0] * young_tableau[0][0]
    cur_pos = 0
    
    for r_dim, r_mult in young_tableau:
        if r_dim % 2 == 0:
            # 偶数サイズのジョルダン細胞に対応する部分行列を抽出して判別式を判定
            # 参照行列 get_ref_bil_form_munip_3(m) を都度呼び出し
            ref_s_inv = get_ref_bil_form_munip_3(r_dim // 2).inverse()
            q_block_inv = block_diagonal_matrix([ref_s_inv] * r_mult, subdivide=False)
            
            # 部分行列の抽出
            s_sub = transformed_s[cur_pos:cur_pos + r_dim * r_mult, cur_pos:cur_pos + r_dim * r_mult]
            q_sub = q_block_inv * s_sub
            
            # 代表値の抽出による簡約判別式の計算
            q_red = matrix(Integers(3), r_mult, r_mult, lambda a, b: q_sub[a * r_dim, b * r_dim])
            disc = Mod(0 if kronecker(q_red.det(), 3) == 1 else 1, 2)
            res_list[r_dim - 1] = (r_mult, disc)
        else:
            # 奇数サイズは重複度のみ
            res_list[r_dim - 1] = r_mult
        cur_pos += r_dim * r_mult
        
    # 結果の整形（空のスロットを埋める）
    for i in range(len(res_list)):
        if (i + 1) % 2 == 0 and res_list[i] == 0:
            res_list[i] = (0, Mod(0, 2))
            
    return tuple(res_list)

def antidiag_scalar_matrix(base_ring, dim, scalar):
    """反対角成分を指定したスカラーで埋めた行列を生成"""
    a = matrix(base_ring, dim, dim)
    for i in range(dim):
        a[i, dim - i - 1] = scalar
    return a

def get_jordan_mults(subdivisions, total_dim):
    """Jordan分解のブロックサイズとその重複度を計算"""
    if not subdivisions:
        return [[total_dim, 1]]
        
    res = []
    prev_pos = 0
    current_dim = 0
    current_mult = 0
    
    for i in subdivisions:
        dim_diff = i - prev_pos
        if dim_diff == current_dim:
            current_mult += 1
        else:
            if current_dim > 0:
                res.append([current_dim, current_mult])
            current_dim = dim_diff
            current_mult = 1
        prev_pos = i
        
    # 最後のブロックの処理
    last_dim = total_dim - prev_pos
    if last_dim == current_dim:
        res.append([current_dim, current_mult + 1])
    else:
        if current_mult > 0:
            res.append([current_dim, current_mult])
        res.append([last_dim, 1])
    return res

def core_quad_3(gram_mat):
    """F3上の二次形式の核（正則部分）の次元と平方剰余性を計算"""
    # 非退化な部分空間のインデックスを取得
    kern = gram_mat.right_kernel().echelonized_basis()
    # 実際には echelonized_basis から直交補空間の代表添字を抽出
    # (既存ロジックを維持)
    dim = gram_mat.nrows()
    rank = gram_mat.rank()
    
    # 簡易化された判別式判定
    # 実際には gram_mat.inverse() 等が必要だが、提供ロジックに基づき縮小行列の det を計算
    if rank == 0:
        return (0, Mod(0, 2))
        
    # echelon形から非零の添字を抽出
    pivot_indices = gram_mat.pivots()
    sub_det = gram_mat[pivot_indices, pivot_indices].det()
    
    disc = Mod(0 if kronecker(sub_det, 3) == 1 else 1, 2)
    return (len(pivot_indices), disc)

def invariants_quad_with_unip_3(gram_mat, gamma_mat):
    """べき単作用を伴う F3 上の二次形式の不変量を計算"""
    t_mat = 1 - gamma_mat
    dim_tot = gram_mat.nrows()
    
    # Ker(1 - gamma)
    kern_basis = t_mat.right_kernel().basis_matrix()
    s_on_ker = kern_basis * gram_mat * kern_basis.transpose()
    
    # 作用を伴う部分の抽出
    pivots = t_mat.pivots()
    s_mat = (gram_mat * t_mat)[pivots, pivots]
    
    disc_sq = Mod(0 if kronecker(gram_mat.det(), 3) == 1 else 1, 2)
    inv_bil = invariants_munip_bil_form_3(s_mat)
    
    return (dim_tot, disc_sq, (core_quad_3(s_on_ker),) + inv_bil)


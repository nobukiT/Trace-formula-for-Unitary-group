class QuadMod4:
    """Unimodular lattice mod 4"""
    def __init__(self, dim, evenpart, oddpart):
        if dim < 0:
            raise ValueError("QuadMod4: negative dimension.")
        self.dim = dim
        self.evenpart = Mod(evenpart, 2)
        self.oddpart = oddpart
        
        # oddpart に基づく奇数部分の次元決定
        if self.oddpart == 0:
            self.odddim = 0
        elif self.oddpart in (1, 2):
            self.odddim = 1
        elif self.oddpart in (3, 4, 5):
            self.odddim = 2
        else:
            raise ValueError(f"QuadMod4: invalid oddpart {oddpart}")
            
        if self.odddim > dim or (self.evenpart != 0 and dim < self.odddim + 2):
            raise ValueError("QuadMod4: inconsistent dimensions.")

    def __neg__(self):
        new_oddpart = 0
        if self.oddpart in (1, 2):
            new_oddpart = 3 - self.oddpart
        elif self.oddpart in (3, 4, 5):
            new_oddpart = 8 - self.oddpart
        return QuadMod4(self.dim, self.evenpart, new_oddpart)

    def __eq__(self, other):
        if not isinstance(other, QuadMod4): return False
        return (self.dim == other.dim and 
                self.evenpart == other.evenpart and 
                self.oddpart == other.oddpart)

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        res = f"Unimodular lattice mod 4, dimension {self.dim}"
        if self.dim > 0:
            parts = []
            if self.evenpart == 0 and self.dim > self.odddim:
                parts.append(f"{(self.dim - self.odddim) // 2} hyperbolic planes")
            elif self.evenpart == 1:
                if self.dim > self.odddim + 2:
                    parts.append(f"{(self.dim - self.odddim) // 2 - 1} hyperbolic planes")
                parts.append("[[2,1],[1,2]]")
            
            odd_labels = {
                1: "[1]", 2: "[3]", 3: "[[1,0],[0,1]]", 
                4: "[[1,0],[0,3]]", 5: "[[3,0],[0,3]]"
            }
            if self.oddpart in odd_labels:
                parts.append(odd_labels[self.oddpart])
            
            res += f", (sum of { ' and '.join(parts) })"
        return res

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash((self.dim, self.evenpart, self.oddpart))

def planeorthmod4_inplace(A, i, j):
    """i, j 平面の直交補空間へ基底を射影 (mod 4)"""
    n = A.nrows()
    detinv = A[i, i] * A[j, j] - A[i, j] * A[j, i]
    for k in range(n):
        if k != i and k != j:
            a = detinv * (A[j, j] * A[k, i] - A[i, j] * A[k, j])
            b = detinv * (-A[j, i] * A[k, i] + A[i, i] * A[k, j])
            if a != 0 or b != 0:
                for l in range(n):
                    A[k, l] -= (a * A[i, l] + b * A[j, l])
                    A[l, k] = A[k, l] # 対称性を維持
    return

def QuadMod4_from_mat(S):
    """行列 S から QuadMod4 オブジェクトを抽出"""
    dim = S.nrows()
    if dim != S.ncols():
        raise ValueError("Matrix must be square.")
    
    A = matrix(Integers(4), S)
    n = dim
    evenpart = Mod(0, 2)
    
    # 簡約化ループ
    while n > 2:
        # A[0,0] を偶数にする調整
        if Mod(A[0, 0], 2) != 0:
            if Mod(A[1, 1], 2) != 0:
                A.add_multiple_of_row(0, 1, 1)
                A.add_multiple_of_column(0, 1, 1)
            else:
                A.swap_rows(0, 1)
                A.swap_columns(0, 1)
        
        # 偶基底の探索
        i = 1
        while i < n and (Mod(A[i, i], 2) != 0 or Mod(A[0, i], 2) == 0):
            i += 1
            
        if i == n: # 偶平面が見つからない場合
            i = 1
            while i < n and (Mod(A[i, i], 2) == 0 or Mod(A[0, i], 2) == 0):
                i += 1
            if i == n: raise ValueError("Not unimodular.")
            
            j = 1
            while j < n and (Mod(A[j, j], 2) == 0 or Mod(A[0, j], 2) != 0):
                j += 1
            
            if j == n: # 直交補空間が偶
                A.swap_rows(i, n - 1); A.swap_columns(i, n - 1)
                A.swap_rows(0, n - 2); A.swap_columns(0, n - 2)
                planeorthmod4_inplace(A, n - 2, n - 1)
                # 再探索
                i = 1
                while i < n and (Mod(A[i, i], 2) != 0 or Mod(A[0, i], 2) == 0): i += 1
            else:
                A.add_multiple_of_row(i, j, 1); A.add_multiple_of_column(i, j, 1)
        
        if i > 1:
            A.swap_rows(1, i); A.swap_columns(1, i)
        
        planeorthmod4_inplace(A, 0, 1)
        if A[0, 0] == 2 and A[1, 1] == 2:
            evenpart += 1
        
        n -= 2
        A = A.submatrix(2, 2)

    # 最終的な判定
    if n == 2 and Mod(A[0, 0], 2) == 0 and Mod(A[1, 1], 2) == 0:
        if A[0, 0] == 2 and A[1, 1] == 2:
            evenpart += 1
        n = 0

    if n == 0:
        oddpart = 0
    elif n == 1:
        oddpart = 1 if A[0, 0] == 1 else 2
    else: # n == 2
        if evenpart == 1 and A.det() == 1:
            evenpart = 0
            A = -A
        
        # 1 を表現するか？
        if A[0, 0] == 1 or A[1, 1] == 1 or (A[0, 0] + A[1, 1] + 2 * A[0, 1]) % 4 == 1:
            oddpart = 3 if A.det() == 1 else 4
        else:
            oddpart = 5
            
    return QuadMod4(dim, evenpart, oddpart)

def planeorthmod4_inplace_newbasis(A, n, P, i, j):
    """
    インデックス i, j で指定される平面を、行列 A の他の基底から直交化する。
    A: Z/4Z 上の対称行列 (破壊的変更)
    n: 行列の次数
    P: Z/2Z 上の基底変換行列 (破壊的変更)
    """
    # Z/4Z における平面 (i, j) の行列式 (あるいはその逆元) を計算
    # ユニモジュラ性の仮定から、det は 1 または 3 (=-1) であり、det == det_inv
    det_inv = A[i, i] * A[j, j] - A[i, j] * A[j, i]
    
    # i, j 以外のすべての基底ベクトル k に対して直交化を適用
    for k in range(n):
        if k == i or k == j:
            continue
            
        # 平面 (i, j) の逆行列を利用して、k 番目のベクトルを平面に射影した際の係数を求める
        # a, b は A[k] = a*A[i] + b*A[j] (mod 平面の直交補空間) となる係数
        a = det_inv * (A[j, j] * A[k, i] - A[i, j] * A[k, j])
        b = det_inv * (-A[j, i] * A[k, i] + A[i, i] * A[k, j])
        
        # A の k 列目と k 行目から、i, j の成分を引いて直交させる
        A.add_multiple_of_column(k, i, -a)
        A.add_multiple_of_column(k, j, -b)
        A.add_multiple_of_row(k, i, -a)
        A.add_multiple_of_row(k, j, -b)
        
        # 変換行列 P (Z/2Z) にもこの操作を反映させる
        # a, b を Z/2Z に落として加算
        if Mod(a, 2) != 0:
            P.add_multiple_of_column(k, i, 1)
        if Mod(b, 2) != 0:
            P.add_multiple_of_column(k, j, 1)

def QuadMod4_newbasis_from_mat(S):
    """
    ユニモジュラ格子 S を Z/4Z 上で正規形に変換し、不変量と変換行列を返す。
    """
    dim = S.nrows()
    if dim != S.ncols():
        raise ValueError("Error in quadmod4_newbasis: nonsquare matrix cannot be symmetric.")

    # 計算用の環を設定
    R4 = Integers(4)
    R2 = Integers(2)
    
    A = matrix(R4, S)
    new_basis = identity_matrix(R2, dim)
    n = dim
    even_part = Mod(0, 2)

    # 1. 偶ユニモジュラ平面 (Even unimodular planes) の抽出ループ
    while n > 2 or (n == 2 and Mod(A[0, 0], 2) == 0 and Mod(A[1, 1], 2) == 0):
        # 最後の基底ベクトル A[n-1] を偶(even)にする
        if Mod(A[n-1, n-1], 2) != 0:
            if Mod(A[n-2, n-2], 2) != 0:
                # e_n を e_n + e_{n-1} で置き換え
                A[n-1, n-1] += 2 * A[n-1, n-2] + A[n-2, n-2]
                for i in range(n - 1):
                    val = A[n-2, i]
                    A[n-1, i] += val
                    A[i, n-1] += val
                new_basis.add_multiple_of_column(n-1, n-2, 1)
            else:
                # e_{n-1} と e_n を入れ替え
                A.swap_columns(n-2, n-1)
                A.swap_rows(n-2, n-1)
                new_basis.swap_columns(n-2, n-1)

        # e_n と直交しない偶基底ベクトルを探す
        idx = -1
        for i in range(n - 2, -1, -1):
            if Mod(A[i, i], 2) == 0 and Mod(A[n-1, i], 2) != 0:
                idx = i
                break
        
        if idx < 0:
            # 偶基底が見つからない場合、奇基底を探す
            for i in range(n - 2, -1, -1):
                if Mod(A[i, i], 2) != 0 and Mod(A[n-1, i], 2) != 0:
                    idx = i
                    break
            
            if idx < 0:
                raise RuntimeError("Error: quadmod4_newbasis: not unimodular.")

            # e_n と直交する奇基底ベクトルを探す
            j_idx = -1
            for j in range(n - 2, -1, -1):
                if Mod(A[j, j], 2) != 0 and Mod(A[n-1, j], 2) == 0:
                    j_idx = j
                    break
            
            if j_idx < 0:
                # (e_i, e_n) が奇ユニモジュラ平面。直交補空間を偶にする処理。
                if idx > 0:
                    A.swap_columns(idx, 0); A.swap_rows(idx, 0)
                    new_basis.swap_columns(idx, 0)
                A.swap_columns(1, n-1); A.swap_rows(1, n-1)
                new_basis.swap_columns(1, n-1)
                
                # 平面直交化（外部関数を想定）
                planeorthmod4_inplace_newbasis(A, n, new_basis, 0, 1)
                
                # 再度検索
                for i in range(n - 2, -1, -1):
                    if Mod(A[i, i], 2) == 0 and Mod(A[n-1, i], 2) != 0:
                        idx = i
                        break
                if idx < 0:
                    raise RuntimeError("Error: quadmod4_newbasis reduction failed.")
            else:
                # (e_i + e_j, e_n) は偶ユニモジュラ平面になる
                A.add_multiple_of_column(idx, j_idx, 1)
                A.add_multiple_of_row(idx, j_idx, 1)
                new_basis.add_multiple_of_column(idx, j_idx, 1)

        # (e_idx, e_{n-1}) を偶ユニモジュラ平面として確定させる
        if idx < n - 2:
            A.swap_columns(n-2, idx); A.swap_rows(n-2, idx)
            new_basis.swap_columns(n-2, idx)
        
        planeorthmod4_inplace_newbasis(A, n, new_basis, n-2, n-1)

        # 標準形 diag(0,2) または diag(2,0) の調整
        if A[n-2, n-2] == 0 and A[n-1, n-1] == 2:
            A.add_multiple_of_column(n-1, n-2, 1); A.add_multiple_of_row(n-1, n-2, 1)
            new_basis.add_multiple_of_column(n-1, n-2, 1)
        elif A[n-2, n-2] == 2 and A[n-1, n-1] == 0:
            A.add_multiple_of_column(n-2, n-1, 1); A.add_multiple_of_row(n-2, n-1, 1)
            new_basis.add_multiple_of_column(n-2, n-1, 1)

        # even_part の更新
        if even_part == 0:
            if A[n-1, n-1] == 2:
                even_part = Mod(1, 2)
        else:
            if A[n-1, n-1] == 0:
                # 平面の入れ替え
                indices = [n-2, n-1, n, n+1]
                for k in range(0, 2):
                    A.swap_rows(indices[k], indices[k+2])
                    A.swap_columns(indices[k], indices[k+2])
                    new_basis.swap_columns(indices[k], indices[k+2])
            else:
                even_part = Mod(0, 2)
                # 基底のクリアと更新
                for k in [n-2, n-1, n, n+1]: A[k, k] = 0
                new_basis.add_multiple_of_column(n-2, n, 1)
                new_basis.add_multiple_of_column(n+1, n-1, 1)
                new_basis.add_multiple_of_column(n-1, n-2, 1)
                new_basis.add_multiple_of_column(n, n+1, 1)
        
        n -= 2

    # 2. 残った奇部分 (Odd part) の分類
    odd_part = 0
    if n == 1:
        odd_part = 1 if A[0, 0] == 1 else 2
    elif n == 2:
        # 対角成分が 0 の場合の調整
        if Mod(A[0, 0], 2) == 0:
            A.add_multiple_of_column(0, 1, 1); A.add_multiple_of_row(0, 1, 1)
            new_basis.add_multiple_of_column(0, 1, 1)
        elif Mod(A[1, 1], 2) == 0:
            A.add_multiple_of_column(1, 0, 1); A.add_multiple_of_row(1, 0, 1)
            new_basis.add_multiple_of_column(1, 0, 1)

        # 特殊ケースの結合
        if even_part == 1 and A[0, 0] * A[1, 1] == 1:
            even_part = 0
            # A[0,0], A[1,1] と A[2,2], A[3,3] の相互作用を解消する操作
            for i in [0, 1]:
                for j in [2, 3]:
                    A.add_multiple_of_column(i, j, 1)
                    A.add_multiple_of_column(j, i, 1)
                    A.add_multiple_of_row(i, j, 1)
                    A.add_multiple_of_row(j, i, 1)
                    new_basis.add_multiple_of_column(i, j, 1)
                    new_basis.add_multiple_of_column(j, i, 1)

        # 奇部分の不変量を決定
        if A[0, 0] == 3:
            if A[1, 1] == 1:
                new_basis.swap_columns(0, 1)
                odd_part = 4 # diag(1,3)
            else:
                odd_part = 5 # diag(3,3)
        else:
            odd_part = 3 if A[1, 1] == 1 else 4 # diag(1,1) or diag(1,3)

    return QuadMod4(dim, even_part, odd_part), new_basis

class QuadQp:
    """Quadratic form over Q_p"""
    def __init__(self, p, dim, disc, hasse):
        if not is_prime(p): raise ValueError(f"{p} is not prime")
        if dim <= 0: raise ValueError("Dimension must be >= 1")
        if disc == 0 or hasse**2 != 1: raise ValueError("Invalid disc or Hasse invariant")
        
        self.p = p
        self.dim = dim
        self.disc = self._unique_repr_mod_squares(disc, p)
        self.hasse = hasse
        
        # 次元の整合性チェック
        if dim == 1 and hasse == -1:
            raise ValueError("Impossible form in dimension 1")
        if dim == 2 and hasse == -1:
            if (p > 2 and kronecker(-self.disc, p) == 1) or (p == 2 and Mod(self.disc, 8) == 7):
                raise ValueError("Impossible form in dimension 2")

    @staticmethod
    def _unique_repr_mod_squares(d, p):
        v = d.valuation(p)
        a = d // p**v
        if p == 2:
            for i in (1, 3, 5, 7):
                if Mod(a - i, 8) == 0: break
        else:
            for i in range(1, p):
                if kronecker(a * i, p) == 1: break
        return i * p if v % 2 == 1 else i

    def has_disc(self, d):
        v = d.valuation(self.p)
        if (v - self.disc.valuation(self.p)) % 2 != 0:
            return False
        d_unit = d // self.p**v if v % 2 == 0 else d // self.p**(v - 1)
        if self.p > 2:
            return kronecker(self.disc / d_unit, self.p) == 1
        else:
            return Mod(self.disc / d_unit, 8) == 1

    def __add__(self, other):
        if self.p != other.p: raise ValueError("Primes do not match.")
        new_hasse = self.hasse * other.hasse * hilbert_symbol(self.disc, other.disc, self.p)
        return QuadQp(self.p, self.dim + other.dim, self.disc * other.disc, new_hasse)

    def contains_regular_lattice(self):
        if self.dim == 0: return True
        if self.p != 2:
            return Mod(self.disc.valuation(self.p), 2) == 0 and self.hasse == 1
        else:
            n = self.dim // 2
            if self.dim % 2 == 1:
                v = self.disc.valuation(2)
                if v % 2 != 1: return False
                x = ((-1)**n * self.disc) // 2**v
                return self.hasse == (-1)**(n * (n - 1) // 2) * hilbert_symbol(x, -1, 2)**n
            else:
                if self.has_disc((-1)**n):
                    return self.hasse == (-1)**(n * (n - 1) // 2)
                elif self.has_disc(3 * (-1)**(n - 1)):
                    return self.hasse == (-1)**(1 + n * (n - 1) // 2)
        return False

    def __str__(self):
        return (f"Quadratic form over QQ_{self.p} in dim {self.dim} "
                f"with disc {self.disc} and Hasse {self.hasse}")

    def __repr__(self):
        return str(self)

def quadQp_from_mat(A, p):
    dim = A.nrows()
    
    # 行列のベースリングが QQ でない（Qpなど）場合、有理数にリフトする
    if A.base_ring() != QQ:
        # p進数の要素は .lift() で有理数(ZZ)の表現に戻せる
        A_QQ = matrix(QQ, dim, dim, [x.lift() if hasattr(x, 'lift') else QQ(x) for x in A.list()])
    else:
        A_QQ = A
        
    disc = A_QQ.det()
    
    if disc == 0:
        raise ValueError(f"Discriminant is zero at prime {p}. Matrix:\n{A_QQ}")
    
    # QuadraticForm は必ず QQ 上で宣言する (hilbert_symbol エラー回避のため)
    Q = QuadraticForm(QQ, 2 * A_QQ)
    hasse = Q.hasse_invariant(p)
    
    # hasse不変量を厳格に 1 か -1 にキャスト
    return QuadQp(p, dim, disc, 1 if hasse > 0 else -1)

def sum_list_quad(L):
    """QuadQp リストの和"""
    if not L: raise ValueError("Empty list.")
    res = L[0]
    for item in L[1:]:
        res += item
    return res

class QuadGlob:
    """Global Quadratic form over QQ"""
    def __init__(self, dim, disc, conductor, negdim):
        self.dim = dim
        self.disc = self._square_reduce(disc)
        self.conductor = self._square_reduce(conductor)
        self.negdim = negdim
        
        if dim <= 0 or negdim < 0 or negdim > dim:
            raise ValueError("Invalid QuadGlob parameters.")
        if ((-1)**negdim * self.disc < 0 or 
            (-1)**(negdim * (negdim - 1) // 2) * moebius(self.conductor) != 1):
            raise ValueError("Inconsistent QuadGlob data.")

    @staticmethod
    def _square_reduce(n):
        if n == 0: return 0
        if n < 0: return -QuadGlob._square_reduce(-n)
        res = n
        d = 2
        while d * d <= res:
            while res % (d * d) == 0:
                res //= (d * d)
            d += 1
        return res

    def __str__(self):
        return (f"Quadratic form over QQ: signature ({self.dim - self.negdim}, {self.negdim}), "
                f"disc {self.disc}, conductor {self.conductor}")

    def __repr__(self):
        return str(self)
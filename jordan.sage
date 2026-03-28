def jordanhermi_unr_generalized(A, E, OE, residue_field, unif, conj, val_func):
    """
    Jordan decomposition of a Hermitian matrix over a general local algebraic number field (unramified case).
    
    Arguments:
      A: The Hermitian matrix to be Jordan decomposed
      E: The field to which the components belong (quadratic extension E/F)
      OE: The ring of integers of E
      residue_field: The residue field of E (F_q)
      unif: A prime element of the base field F (uniformizer varpi_F)
      conj: A non-trivial automorphism (conjugation) of the Galois group of E/F
      val_func: A function that takes an element of E and returns the valuation at the corresponding prime ideal
    """
    dim = A.nrows()
    d = 0
    dims = []
    B = copy(A)
    newbasis = Matrix(E, dim, dim, 1)

    # p ではなく、与えられた付値関数 val_func を使用して最小付値を計算
    vals_B = [val_func(x) for x in B.list() if x != 0]
    i = min(vals_B) if vals_B else 0

    while d < dim:
        # p**i ではなく基礎体の素元 unif**i でスケール
        while Matrix(residue_field, B / (unif**i)) == 0:
            i += 1
            dims.append(0)

        Bmod = Matrix(residue_field, B / (unif**i))
        rank = Bmod.rank()
        dims.append(rank)

        kernbase = Bmod.right_kernel().echelonized_basis()
        indices = []
        k = 0
        for v in kernbase:
            l = k + 1
            while v[l-1] == 0:
                l += 1
            indices.extend(range(k+1, l))
            k = l
        indices.extend(range(k+1, dim - d + 1))

        P = Matrix(OE, dim, dim)
        k = 1
        for j in indices:
            P[j-1, k-1] = 1
            k += 1

        # 注意: integer_kernel は SageMath の整数環・代数体サポートに依存します
        kernbase_oe = (B * P).integer_kernel(OE).basis()
        for v in kernbase_oe:
            l = 0
            for x in v:
                P[l, k-1] = conj(x)
                l += 1
            k += 1

        if d + rank < dim:
            Ptransconj = Matrix(OE, [[conj(P[b, a]) for b in range(dim)] for a in range(dim)])
            B = (Ptransconj * B * P).submatrix(rank, rank)

        newbasis = newbasis * block_diagonal_matrix(diagonal_matrix([1]*d), P, subdivide=False)
        d += rank
        i += 1

    newbasistransconj = Matrix(E, [[conj(newbasis[b, a]) for b in range(dim)] for a in range(dim)])
    return newbasis, newbasistransconj * A * newbasis, dims


def jordanhermi_split_generalized(A, E, OE, residue_field, unif, conj, val_func, t):
    """
    一般の局所代数体におけるエルミート行列のジョルダン分解（分裂ケース）。
    
    追加引数:
      t: E が 2つの素イデアル P1, P2 に分裂するとき、
         t = 1 mod P1 かつ t = 0 mod P2 を満たす OE の元（局所的な冪等元の近似）
    """
    tconj = OE(conj(t))
    
    dim = A.nrows()
    d = 0
    dims = []
    B = copy(A)
    newbasis = Matrix(E, dim, dim, 1)

    vals_B = [val_func(x) for x in B.list() if x != 0]
    i = min(vals_B) if vals_B else 0

    while d < dim:
        while Matrix(residue_field, B / (unif**i)) == 0:
            i += 1
            dims.append(0)

        Bmod = Matrix(residue_field, B / (unif**i))
        rank = Bmod.rank()
        dims.append(rank)

        rkern = Bmod.right_kernel().echelonized_basis()
        indices = []
        k = 0
        for v in rkern:
            l = k + 1
            while v[l-1] == 0:
                l += 1
            indices.extend(range(k+1, l))
            k = l
        indices.extend(range(k+1, dim - d + 1))
        
        P1 = Matrix(OE, dim, dim)
        k = 1
        for j in indices:
            P1[j-1, k-1] = 1
            k += 1

        lkern = Bmod.left_kernel().echelonized_basis()
        indices = []
        k = 0
        for v in lkern:
            l = k + 1
            while v[l-1] == 0:
                l += 1
            indices.extend(range(k+1, l))
            k = l
        indices.extend(range(k+1, dim - d + 1))
        
        P2 = Matrix(OE, dim, dim)
        k = 1
        for j in indices:
            P2[j-1, k-1] = 1
            k += 1

        # 冪等元 t を用いた基底の局所的な貼り合わせ
        P = (1-t)*P1 + t*P2
        Ptransconj = (1-tconj)*P1.transpose() + tconj*P2.transpose()

        kernbase_oe = (B * P).integer_kernel(OE).basis()
        for v in kernbase_oe:
            l = 0
            for x in v:
                P[l, k-1] = conj(x)
                Ptransconj[k-1, l] = x
                l += 1
            k += 1

        if d + rank < dim:
            B = (Ptransconj * B * P).submatrix(rank, rank)

        newbasis = newbasis * block_diagonal_matrix(diagonal_matrix([1]*d), P, subdivide=False)
        d += rank
        i += 1

    newbasistransconj = Matrix(E, [[conj(newbasis[b, a]) for b in range(dim)] for a in range(dim)])
    return newbasis, newbasistransconj * A * newbasis, dims


def jordanhermi_generalized(A, is_split, E, OE, residue_field, unif, conj, val_func, t=None):
    """
    一般化されたジョルダン分解の統合ラッパー関数
    """
    if is_split:
        if t is None:
            raise ValueError("Split case requires the idempotent approximator 't'.")
        return jordanhermi_split_generalized(A, E, OE, residue_field, unif, conj, val_func, t)
    else:
        return jordanhermi_unr_generalized(A, E, OE, residue_field, unif, conj, val_func)
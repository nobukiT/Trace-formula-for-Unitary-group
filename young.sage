def young_list(M):
    """
    べき零行列 M のジョルダン細胞のサイズ（ヤング図形）を計算する。
    
    数学的背景:
    L[k] = rank(M^k) - rank(M^{k+1}) は、ジョルダン標準形における
    サイズが k+1 以上のブロックの個数に対応する（ヤング図形の「列」の長さ）。
    """
    if not M.is_square():
        raise ValueError("行列は正方行列である必要があります。")
        
    total_dim = M.nrows()
    if total_dim == 0:
        return []

    partition = []
    current_nullity = total_dim - M.rank()  # dim(ker(M))
    prev_nullity = 0
    
    # 行列の累乗を保持する変数
    M_pow = M
    
    # 核の次元の増分がなくなるまでループ
    while current_nullity > prev_nullity:
        # この増分がヤング図形の「列」の長さを構成する
        partition.append(current_nullity - prev_nullity)
        
        # 次の累乗の核の次元を計算
        prev_nullity = current_nullity
        M_pow *= M
        current_nullity = total_dim - M_pow.rank()

    return partition
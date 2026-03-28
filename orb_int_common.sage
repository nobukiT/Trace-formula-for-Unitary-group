def degrees_cyclo(n, primes_list):
    """
    指定された素数リストから作られる整数 m のうち、
    オイラーのφ関数の値 phi(m) が n 以下になるような [m, phi(m)] のリストを生成する。
    """
    res = []
    
    def generate_combinations(prime_idx, current_m, current_phi):
        # ベースケース: すべての素数の指数を決定したらリストに追加
        if prime_idx == len(primes_list):
            res.append([current_m, current_phi])
            return
            
        p = primes_list[prime_idx]
        
        # パターン1: この素数 p を1回も使わない場合 (指数 e = 0)
        generate_combinations(prime_idx + 1, current_m, current_phi)
        
        # パターン2: この素数 p を1回以上使う場合 (指数 e >= 1)
        # p^e を掛けたときの phi の増加分は (p - 1) * p^(e - 1)
        m_multiplier = p
        phi_multiplier = p - 1
        
        # phi(m) が n を超えない範囲で、p を何回も掛けていく
        while current_phi * phi_multiplier <= n:
            generate_combinations(
                prime_idx + 1, 
                current_m * m_multiplier, 
                current_phi * phi_multiplier
            )
            # さらに p をもう1回掛ける (指数を +1 する)
            m_multiplier *= p
            phi_multiplier *= p

    # 初期値: m = 1, phi(1) = 1 から探索スタート
    generate_combinations(0, 1, 1)
    
    # 元のコードと互換性を持たせるため、m の昇順にソートして返す
    return sorted(res, key=lambda x: x[0])

def totally_ramified_conj_classes_unitary(N, degrees=None):
    """
    局所体 Q_p 上の完全分岐を含む共役類（円分多項式の積）を生成する（ユニタリ群用）。
    N: 空間の全次元 (U(n,n) なら N = 2n)
    """
    res = []

    primes_list = list(primes(N + 2))
    
    if degrees is None:
        degrees = degrees_cyclo(N, primes_list)
        
    degrees_rev = list(reversed(degrees))
    primes_rev = reversed(primes_list)
    
    for p in primes_rev:
        min_deg_mult = 2 if p == 2 else p - 1
        valid_degrees = [
            (m, phi_m) for m, phi_m in degrees_rev
            if (phi_m * min_deg_mult <= N) 
            and (m % p != 0) 
            and (m % 4 != 2)
        ]
        
        for m, phi_m in valid_degrees:
            # comp_deg[k] = phi(m * p^k) の値
            comp_deg = [phi_m]
            deg_k = phi_m * (p - 1)
            while deg_k <= N:
                comp_deg.append(deg_k)
                deg_k *= p
                
            max_k = len(comp_deg) - 1

            def generate_partitions(target, k_idx=max_k, current_mults=None):
                if current_mults is None:
                    current_mults = [0] * (max_k + 1)
                    
                if k_idx < 0:
                    # 【修正1】current_mults[1] >= 1 の制限を削除し、完全に k=0 だけの不分岐クラスも許容する
                    if target == 0: 
                        yield tuple((k, current_mults[k]) for k in range(max_k + 1) if current_mults[k] > 0)
                    return
                
                coin_val = comp_deg[k_idx]
                # 【修正2】k_idx == 1 でも最低 0 個からスタートできるようにする
                min_count = 0  
                max_count = target // coin_val
                
                for count in range(min_count, max_count + 1):
                    current_mults[k_idx] = count
                    yield from generate_partitions(target - count * coin_val, k_idx - 1, current_mults)
            
            for dims_tup in generate_partitions(N):
                res.append((p, m, dims_tup))
                
    return res
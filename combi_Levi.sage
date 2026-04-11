# ======================================================================
# Combinatorics for Unitary Groups U(p, q)
# ======================================================================

def enum_Levis_Upq(p_sig, q_sig):
    """
    Enumerates Levi subgroups for U(p, q) that support discrete series.
    Over R, GL_k(C) has discrete series iff k=1.
    The Levi subgroups are of the form GL_1(C)^a x U(p-a, q-a).
    
    Returns: list of [a, p-a, q-a] where a is the number of GL_1(C) factors.
    """
    res = []
    # GL_1(C) requires one positive and one negative dimension, 
    # so 'a' cannot exceed min(p_sig, q_sig).
    for a in range(min(p_sig, q_sig) + 1):
        res.append([a, p_sig - a, q_sig - a])
    return res

def set_partitions_pairs(n, a, whole_set=None):
    """
    Selects 'a' disjoint pairs from a set of 'n' elements.
    Used for the GL_1(C) components in the Weyl group of U(p, q).
    """
    if not whole_set:
        whole_set = list(range(1, n + 1))
        
    if a == 0:
        return [[[], whole_set, Mod(0, 2)]]
        
    if len(whole_set) < 2 * a:
        return []
        
    res = []
    first = whole_set[0]
    
    for i in range(1, len(whole_set)):
        second = whole_set[i]
        pair = [first, second]
        rem_set = whole_set[1:i] + whole_set[i+1:]
        
        sign_shift = i - 1 
        
        for pa in set_partitions_pairs(n - 2, a - 1, rem_set):
            res.append([pair + pa[0], pa[1], pa[2] + Mod(sign_shift, 2)])
            
    return res

def index_from_dom_weight_U(n, dom_weight):
    res = sum([binomial(dom_weight[n - i - 1] + i, i + 1) for i in range(n)])
    return res

def contrib_Levi_Upq(L, ell_tr_Upq, ell_tr_GL1C, abs_max_weight):
    """
    Computes the contribution of a specific Levi subgroup GL_1(C)^a x U(p-a, q-a).
    """
    a = L[0]
    p_sub = L[1]
    q_sub = L[2]
    d = p_sub + q_sub  # Dimension of the remaining unitary block
    n = 2 * a + d      # Total dimension p + q
    
    sigma_repr_list = set_partitions_pairs(n, a)
    res = []
    dom_weight = [0] * n
    
    while dom_weight[0] + n - 1 <= abs_max_weight:
        weight_res = 0
        
        for sigma in sigma_repr_list:
            GL1C_contrib = 1
            
            # Contribution from GL_1(C) pairs
            for i in range(a):
                s1 = sigma[0][2 * i]
                s2 = sigma[0][2 * i + 1]
                
                k1 = dom_weight[s1 - 1] - dom_weight[s2 - 1] - s1 + s2 - 1
                
                if Mod(k1, 2) == Mod(0, 2):
                    GL1C_contrib *= ell_tr_GL1C[k1 // 2]
                else:
                    GL1C_contrib = 0
                    break
                    
            if GL1C_contrib == 0:
                continue
                
            # Contribution from U(p-a, q-a)
            U_d_weight = [dom_weight[sigma[1][i] - 1] + n - d + i + 1 - sigma[1][i] for i in range(d)]
            sign_factor = (-1) ** int(sigma[2])
            
            # ell_tr_Upq must be indexed by the signature or the dimension of the remaining block
            weight_res += sign_factor * GL1C_contrib * ell_tr_Upq[d][index_from_dom_weight_U(d, U_d_weight)]
            
        res.append(weight_res)
        
        dom_weight[n - 1] += 1
        i = n - 1
        while i > 0 and dom_weight[i] > dom_weight[i - 1]:
            dom_weight[i] = 0
            dom_weight[i - 1] += 1
            i -= 1
            
    return res

def EPchar_Upq(p_sig, q_sig, ell_tr_Upq, ell_tr_GL1C, abs_max_weight):
    n = p_sig + q_sig
    num_weights = binomial(abs_max_weight + 1, n) 
    res = [0] * num_weights
    
    for L in enum_Levis_Upq(p_sig, q_sig):
        L_contrib = contrib_Levi_Upq(L, ell_tr_Upq, ell_tr_GL1C, abs_max_weight)
        for i in range(num_weights):
            res[i] += L_contrib[i]
            
    return res
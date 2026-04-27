load("quadclass.sage")
load("hermi_density.sage")
load("young.sage")
load("deep_extensions_quad_3.sage")

def quadlattice_oddp_from_hermi_type(p, k, d, mult, I, precomp_trace_matrices, alg_data):
    K = Qp(p)
    if k > 0:
        comp_mat = companion_matrix(cyclotomic_polynomial(p**k), format='right')
        deg_phi = comp_mat.nrows()
    else:
        deg_phi = 1
        
    base_diag_gram_mat = []
    disc_isnorm = Mod(0, 2)
  
    for idx, L_i in enumerate(I):
        i = idx - d
        n_i = L_i[0]
        disc_type = L_i[1] if len(L_i) > 1 else Mod(0, 2)

        if i % 2 == 0:
            if n_i > 0:
                if disc_type == Mod(0, 2):
                    base_diag_gram_mat.extend([precomp_trace_matrices[i+d][0]] * n_i)
                else:
                    base_diag_gram_mat.extend([precomp_trace_matrices[i+d][0]] * (n_i - 1))
                    base_diag_gram_mat.append(precomp_trace_matrices[i+d][1])
                    disc_isnorm += Mod(1, 2)
        else:
            if n_i > 0:
                base_diag_gram_mat.extend([precomp_trace_matrices[i+d]] * (n_i // 2))
                if (n_i % 4 != 0) and (p % 4 != 1):
                    disc_isnorm += Mod(1, 2)

    diag_gram_mat = []
    for M in base_diag_gram_mat:
        diag_gram_mat.extend([matrix(K, M)] * deg_phi)

    gram_mat0 = block_diagonal_matrix(diag_gram_mat) if diag_gram_mat else matrix(K, 0, 0)
    total_dim = gram_mat0.nrows()
    
    if k > 0:
        mult_adj = total_dim // deg_phi
        if mult_adj > 0:
            gamma_mat0 = block_diagonal_matrix([matrix(K, comp_mat)] * mult_adj, subdivide=False)
        else:
            gamma_mat0 = matrix(K, 0, 0)
    else:
        gamma_mat0 = identity_matrix(K, total_dim)

    P, gram_mat, dims_list = jordanquadZp(gram_mat0, p)
    P_K = matrix(K, P)
    
    if total_dim > 0:
        gamma_mat = P_K.inverse() * gamma_mat0 * P_K
        # gram_mat may not be precisely in QQ, but quadQp_from_mat expects something coercible.
        gram_mat_QQ = matrix(QQ, [[QQ(x.lift()) if hasattr(x, 'lift') else QQ(x) for x in row] for row in gram_mat])
        Qpclass = quadQp_from_mat(gram_mat_QQ, p)
    else:
        gamma_mat = matrix(K, 0, 0)
        Qpclass = QuadQp(p, 0, 1, 1) 

    type_list = []
    cur_pos = 0
    for j, r in enumerate(dims_list):
        if r == 0:
            type_list.append([0, Mod(0, 2)])
            continue
        block = gram_mat[cur_pos:cur_pos+r, cur_pos:cur_pos+r]
        det_val = (block / (p**j)).det()
        if det_val == 0:
            raise RuntimeError("Degenerate matrix encountered in Jordan decomposition.")
        u = det_val.unit_part().residue()
        disc_issquare = kronecker(u.lift(), p)
        type_list.append([r, Mod(0 if disc_issquare == 1 else 1, 2)])
        cur_pos += r

    formatted_I = []
    for L_i in I:
        n_i = int(L_i[0]) if len(L_i) > 0 else 0
        t_i = int(L_i[1]) if len(L_i) > 1 else 0
        formatted_I.append([n_i, t_i])

    q = alg_data['q']
    mass_local_term = hermi_ramoddp_mass_local_term(q, d, mult, formatted_I)
    
    result = quadlattice_oddp(p, total_dim, type_list, Qpclass, gram_mat, gamma_mat, mass_local_term)
    result.hermi_type = I
    return [disc_isnorm, result]

class quadlattice_oddp:
    def __init__(self, p, dim, type_list, Qpclass, gram_mat, gamma_mat, mass_local_term):
        self.p = p
        self.dim = dim
        self.type_list = type_list
        self.Qpclass = Qpclass
        self.gram_mat = gram_mat
        self.gamma_mat = gamma_mat
        self.mass_local_term = mass_local_term

        def modp_matrix(M, p):
            Fp = Integers(p)
            n = M.nrows()
            m = M.ncols()
            out = matrix(Fp, n, m)
            for i in range(n):
                for j in range(m):
                    x = M[i, j]
                    if x == 0: continue
                    v = x.valuation()
                    unit = x * (p**(-v))
                    out[i, j] = unit.residue()
            return out

        if len(type_list) <= 2:
            self.is_simple = True
            self.reduction_dim = type_list[1][0] if len(type_list) > 1 else 0
            self.reduction_disc = type_list[1][1] if len(type_list) > 1 else Mod(0, 2)
            u = type_list[0][0] if len(type_list) > 0 else 0
            red_gamma = modp_matrix(gamma_mat[u:, u:], p)
            self.reduction_young = young_list(red_gamma - 1)
            if p == 3:
                red_gram_block = gram_mat[u:, u:] / 3
                red_gram = modp_matrix(red_gram_block, 3)
                self.reduction_invariants = invariants_quad_with_unip_3(red_gram, red_gamma)
        else:
            self.is_simple = False

def quadlattice_2_from_hermi_type(m, D_K, mult, I, precomp_trace_matrices, alg_data):
    k = Integer(m).valuation(2)
    d = alg_data['d']
    Fi = alg_data['F']
    Ei = alg_data['E']
    embed_F_to_E = alg_data['embed_F_to_E']
    p_ideal = alg_data['p_ideal']

    if k > 0:
        comp_mat = companion_matrix(cyclotomic_polynomial(2**k), format='right')
        deg_phi = comp_mat.nrows()
    else:
        deg_phi = 1

    base_diag_gram_mat = []
    disc_is_norm = Mod(0, 2)
    
    for idx, L_i in enumerate(I):
        i = idx - d
        n_i = L_i[0]
        
        if n_i > 0:
            if i % 2 == 1: 
                if k == 2 and n_i % 2 == 1:
                    disc_is_norm += 1
                if len(L_i) > 1 and L_i[1] != 0:
                    base_diag_gram_mat.extend([precomp_trace_matrices[i + d][0]] * (n_i - 1))
                    base_diag_gram_mat.append(precomp_trace_matrices[i + d][1])
                    disc_is_norm += 1
                else:
                    base_diag_gram_mat.extend([precomp_trace_matrices[i + d][0]] * n_i)
            else:
                eps = L_i[1] if len(L_i) > 1 else 0
                eta = L_i[2] if len(L_i) > 2 else 0

                if eps == 0: 
                    base_diag_gram_mat.extend([precomp_trace_matrices[i + d][0]] * (n_i // 2))
                    if k == 2 and n_i % 4 == 2:
                        disc_is_norm += 1
                else:
                    base_diag_gram_mat.extend([precomp_trace_matrices[i + d][0]] * (n_i // 2))
                    if k == 2 and (n_i % 4 in [3, 0]):
                        disc_is_norm += 1
                    
                    if n_i % 2 == 0:
                        base_diag_gram_mat.append(precomp_trace_matrices[i + d][1])
                        
                    if eta == 0:
                        base_diag_gram_mat.append(precomp_trace_matrices[i + d][1])
                    else:
                        base_diag_gram_mat.append(precomp_trace_matrices[i + d][2])
                        disc_is_norm += 1

    diag_gram_mat = []
    for M in base_diag_gram_mat:
        diag_gram_mat.extend([M] * deg_phi)
                        
    gram_mat0 = block_diagonal_matrix(diag_gram_mat) if diag_gram_mat else matrix(QQ, 0, 0)
    total_dim = gram_mat0.nrows()
    
    if k > 0:
        mult_adj = total_dim // deg_phi
        if mult_adj > 0:
            gamma_mat0 = block_diagonal_matrix([comp_mat] * mult_adj, subdivide=False)
        else:
            gamma_mat0 = matrix(QQ, 0, 0)
    else:
        gamma_mat0 = identity_matrix(QQ, total_dim)
    
    P, gram_mat, dims_list = jordanquadZp(gram_mat0, 2)
    P_QQ = matrix(QQ, P)
    
    if total_dim > 0:
        gamma_mat = P_QQ.inverse() * matrix(QQ, gamma_mat0) * P_QQ
        Qpclass = quadQp_from_mat(gram_mat, 2)
    else:
        gamma_mat = matrix(QQ, 0, 0)
        Qpclass = QuadQp(2, 0, 1, 1)
        
    type_list = []
    cur_pos = 0
    
    def is_even(M):
        for r in range(M.nrows()):
            if M[r, r].valuation(2) < 1:
                return False
        return True

    for j, d_size in enumerate(dims_list):
        if d_size == 0:
            type_list.append([0, Mod(0, 2)])
        else:
            if is_even(gram_mat[cur_pos : cur_pos + d_size, cur_pos : cur_pos + d_size] / 2**j):
                type_list.append([d_size, Mod(0, 2)])
            else:
                type_list.append([d_size, Mod(1, 2)])
            cur_pos += d_size

    formatted_I = []
    for L_i in I:
        n_i = int(L_i[0]) if len(L_i) > 0 else 0
        eps = int(L_i[1]) if len(L_i) > 1 else 0
        eta = int(L_i[2]) if len(L_i) > 2 else 0
        formatted_I.append([n_i, eps, eta])

    if Fi == QQ:
        q = p_ideal
    else:
        q = p_ideal.norm()

    mass_local_term = hermi_ram2_mass_local_term(q, d, mult, formatted_I)
    
    result = quadlattice_2(total_dim, type_list, Qpclass, gram_mat, gamma_mat, mass_local_term)
    result.hermi_type = I
    return [disc_is_norm, result]

class quadlattice_2:
    def __init__(self, dim, type_list, Qpclass, gram_mat, gamma_mat, mass_local_term):
        self.dim = dim
        self.type_list = type_list
        self.Qpclass = Qpclass
        self.gram_mat = gram_mat
        self.gamma_mat = gamma_mat
        self.mass_local_term = mass_local_term
        
        if len(type_list) <= 2:
            self.is_simple = True
            if len(type_list) <= 1:
                self.reduction_dim = 0
                self.reduction_young = []
                self.reduction_gamma = matrix(Integers(2), 0, 0)
                self.reduction_mod4 = QuadMod4(0, Mod(0, 2), 0)
                self.reduction_gram = matrix(Integers(4), 0, 0)
            else:
                u = type_list[0][0]
                self.reduction_dim = type_list[1][0]
                S_rat = gram_mat[u:, u:] / 2
                n_sub = S_rat.nrows()
                S = matrix(Integers(4), n_sub, n_sub)
                for r in range(n_sub):
                    for c in range(n_sub):
                        val = S_rat[r, c]
                        try:
                            S[r, c] = Integers(4)(val)
                        except (TypeError, ValueError, ZeroDivisionError):
                            S[r, c] = 0
                g = matrix(Integers(2), gamma_mat[u:, u:])
                self.reduction_mod4, newbasis = QuadMod4_newbasis_from_mat(S)
                try:
                    if g.parent() == newbasis.parent() and newbasis.determinant() % 2 != 0:
                        self.reduction_gamma = newbasis.inverse() * g * newbasis
                    else:
                        self.reduction_gamma = g
                except (ZeroDivisionError, ArithmeticError):
                    self.reduction_gamma = g
                self.reduction_young = young_list(g - 1)
                newbasis_lifted = matrix(QQ, [[QQ(x.lift()) if hasattr(x, 'lift') else QQ(x) for x in row] for row in newbasis])
                S_lifted = matrix(QQ, [[QQ(x.lift()) if hasattr(x, 'lift') else QQ(x) for x in row] for row in S])
                gram_lifted = newbasis_lifted.transpose() * S_lifted * newbasis_lifted
                self.reduction_gram = matrix(Integers(4), gram_lifted)
                if len(self.reduction_young) < 2:
                    self.reduction_hs_mod4 = hyperspacesmod4(S, g)
        else:
            self.is_simple = False

    def __str__(self):
        return (f"Quadratic lattice over ZZ_2, rank {self.dim}, type {self.type_list}, "
                f"disc in QQ_2: {self.Qpclass.disc}, Hasse-Witt: {self.Qpclass.hasse}\n"
                f"and local term in the mass formula: {self.mass_local_term}.")

    def __repr__(self):
        return str(self)

def hyperspacesmod4(S, g):
    n = S.nrows()
    if n == 0: return {}
    if n == 1: return {QuadMod4_from_mat(matrix(Integers(4), 0, 0)): 1}
    res = {}
    X = matrix(Integers(2), 1, n)
    X[0, n-1] = 1
    normX = Mod(S[n-1, n-1], 2)
    i = n - 1
    while i >= 0:
        if normX != 0 and (g * X.transpose() == X.transpose()):
            new_gram = lineorthmod4(S, X)
            Q = QuadMod4_from_mat(new_gram)
            res[Q] = res.get(Q, 0) + 1
        i = n - 1
        while i >= 0 and X[0, i] != 0:
            normX += Mod(S[i, i], 2)
            X[0, i] = 0
            i -= 1
        if i >= 0:
            normX += Mod(S[i, i], 2)
            X[0, i] = 1
    return res

def lineorthmod4(S, X):
    Pmod2 = (X * S).right_kernel().basis_matrix().transpose()
    P = matrix(Integers(4), Pmod2.lift())
    return P.transpose() * S * P
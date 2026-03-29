load("hermi_density.sage")
load("young.sage")

def build_hermilattice_unr(p, k, m2, dim, I, alg_data, t=None):
    """
    Computes the data and local volume of an unramified (or totally ramified) 
    Hermitian lattice over a general local algebraic field.

    INPUT:
        p                      : int; Characteristic of the residue field.
        k                      : int; Exponent for the cyclotomic extension (p^k).
        m2                     : int; Order of the cyclotomic field.
        dim                    : int; Total dimension of the lattice.
        I                      : list; Invariants specifying the Jordan decomposition.
        precomp_trace_matrices : list; Precomputed trace matrices for efficiency.
        alg_data               : dict; Dictionary containing local algebraic structures.
        t                      : (Optional) Parameter passed to the generalized Jordan decomposition.

    OUTPUT:
        list: [disc_isnorm, lattice_obj]
    """
    diag_gram_mat = []
    
    if k > 0:
        deg_phi = cyclotomic_polynomial(p**k).degree()
    else:
        deg_phi = 1
        
    for i, count in enumerate(I):
        if count > 0:
            val = alg_data['unif']**(i - alg_data['d'])
            diag_gram_mat.extend([val] * (count * deg_phi))
            
    hermi_dim_expanded = len(diag_gram_mat)
            
    if k > 0:
        gram_mat0 = diagonal_matrix(alg_data['E'], diag_gram_mat)
        newbasis, gram_mat, type_list = jordanhermi(gram_mat0, False, alg_data, t)
        
        comp_mat_ZZ = companion_matrix(cyclotomic_polynomial(p**k), format='right')
        comp_mat = matrix(alg_data['E'], comp_mat_ZZ)
        
        num_blocks = hermi_dim_expanded // deg_phi
        gamma_mat0 = block_diagonal_matrix([comp_mat] * num_blocks, subdivide=False)
        
        gamma_mat = newbasis**(-1) * gamma_mat0 * newbasis
    else:
        gram_mat = diagonal_matrix(alg_data['E'], diag_gram_mat)
        gamma_mat = identity_matrix(alg_data['E'], hermi_dim_expanded)
        type_list = I

    mass_local_term = hermi_unr_mass_local_term(alg_data['q'], dim, I)
    disc_isnorm = sum((i - alg_data['d']) * count for i, count in enumerate(I)) % 2
        
    lattice_obj = hermilattice_unr(alg_data['q'], p, m2, dim * alg_data['e'], type_list, gram_mat, gamma_mat, mass_local_term)
    
    return [disc_isnorm, lattice_obj]


class hermilattice_unr:
    """
    A data class representing an unramified Hermitian lattice over a general local field.
    """
    def __init__(self, q, p, m2, dim, type_list, gram_mat, gamma_mat, mass_local_term, compute_reduction=True):
        self.q = q  # Order of the residue field
        self.p = p  # Characteristic
        self.m2 = m2
        self.dim = dim
        self.type_list = type_list
        self.gram_mat = gram_mat
        self.gamma_mat = gamma_mat
        self.mass_local_term = mass_local_term
        
        self.is_simple = False

        # 3. Calculate the Young diagram (Jordan normal form) over the residue field
        if compute_reduction and len(type_list) <= 2:
            self.is_simple = True
            
            if len(type_list) <= 1:
                self.reduction_dim = 0
                self.reduction_young = []
            else:
                self.reduction_dim = type_list[1]
                residue_field_local = GF(q, 'a')
                sub_gamma = gamma_mat.submatrix(type_list[0], type_list[0])
                try:
                    reduced_matrix = matrix(residue_field_local, self.reduction_dim, self.reduction_dim, sub_gamma) - 1
                except TypeError:
                    reduced_matrix = matrix(residue_field_local, self.reduction_dim, self.reduction_dim, 
                                            [residue_field_local(x) for x in sub_gamma.list()]) - 1
                self.reduction_young = young_list(reduced_matrix)

    def __str__(self):
        return (f"Hermitian lattice over local field with residue field F_{self.q}, "
                f"rank {self.dim}, type {self.type_list}, "
                f"local term in mass formula: {self.mass_local_term}")

    def __repr__(self):
        return str(self)
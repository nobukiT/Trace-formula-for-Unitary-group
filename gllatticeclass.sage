load("gl_density.sage")
load("young.sage")

def gl_lattice_unr_from_type(p, k, m2, d, e, dim, I):
    
    if k > 0:
        comp_mat = companion_matrix(cyclotomic_polynomial(p**k), format='right')
        gamma_mat = block_diagonal_matrix([comp_mat] * dim, subdivide=False)
        type_list = I
    else:
        gamma_mat = identity_matrix(dim)
        type_list = I

    mass_local_term = gl_mass_local_term(p, k, m2, I)
    
    disc_dummy = Mod(0, 2) 
    
    lattice_obj = gl_lattice_unr(p, m2, dim * e, type_list, gamma_mat, mass_local_term)
    return [disc_dummy, lattice_obj]


class gl_lattice_unr:
    def __init__(self, p, m2, dim, type_list, gamma_mat, mass_local_term, compute_reduction=True):
        self.p = p
        self.m2 = m2
        self.dim = dim
        self.type_list = type_list
        self.GLclass = tuple(type_list)
        self.gamma_mat = gamma_mat
        self.mass_local_term = mass_local_term
        self.is_simple = False
        
        if compute_reduction and len(type_list) <= 2:
            self.is_simple = True
            if len(type_list) <= 1:
                self.reduction_dim = 0
                self.reduction_young = []
            else:
                self.reduction_dim = type_list[1]
                K = CyclotomicField(m2)
                residue_field = K.residue_field(K.primes_above(p)[0])
                sub_gamma = gamma_mat.submatrix(type_list[0], type_list[0])
                reduced_matrix = matrix(residue_field, self.reduction_dim, self.reduction_dim, sub_gamma) - 1
                self.reduction_young = young_list(reduced_matrix)

    def __str__(self):
        return (f"GL lattice over ZZ_{self.p}[zeta_{self.m2}], "
                f"rank {self.dim}, type {self.type_list}, "
                f"local term in mass formula: {self.mass_local_term}")

    def __repr__(self):
        return str(self)
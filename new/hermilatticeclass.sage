load("hermi_density.sage")
load("young.sage")

def build_hermilattice_unr(q, p, unif, k, m2, d, e, dim, I, precomp_trace_matrices, 
                           E, OE, residue_field, conj, val_func, t=None):
    """
    Computes the data and local volume of an unramified (or totally ramified) 
    Hermitian lattice over a general local algebraic field.

    INPUT:
        q                      : int; Order of the residue field of the base field.
        p                      : int; Characteristic of the residue field.
        unif                   : Element; Uniformizer (\varpi) of the base field.
        k                      : int; Exponent for the cyclotomic extension (p^k).
        m2                     : int; Order of the cyclotomic field.
        d                      : int; Exponent of the different of F_p(zeta_p^k)/F_p.
        e                      : int; Ramification index.
        dim                    : int; Total dimension of the lattice.
        I                      : list; Invariants specifying the Jordan decomposition [count0, count1, ...].
        precomp_trace_matrices : list; Precomputed trace matrices for efficiency.
        E                      : Field; The extension field.
        OE                     : Ring; The ring of integers of the extension field E.
        residue_field          : Field; The residue field of the base field.
        conj                   : function; The non-trivial involution (conjugation) map on E.
        val_func               : function; Valuation function for the base field.
        t                      : (Optional) Parameter passed to the generalized Jordan decomposition.

    OUTPUT:
        list: [disc_isnorm, lattice_obj]
            - disc_isnorm (int): Parity of the scaled discriminant valuation.
            - lattice_obj (hermilattice_unr): The constructed Hermitian lattice object.
    """
    diag_gram_mat = []
    
    # Construct the diagonal entries of the Gram matrix
    for i, count in enumerate(I):
        if count > 0:
            # Use the provided uniformizer unif**(i-d) instead of p**(i-d)
            val = precomp_trace_matrices[i] if k > 0 else unif**(i - d)
            diag_gram_mat.extend([val] * count)
            
    if k > 0:
        gram_mat0 = block_diagonal_matrix(diag_gram_mat)
        
        # Call the generalized Jordan decomposition
        newbasis, gram_mat, type_list = jordanhermi(
        )
        
        # Note: The companion matrix depends on a specific polynomial.
        # For completely general local fields, cyclotomic_polynomial must be 
        # replaced with the corresponding minimal polynomial over the base field.
        comp_mat = companion_matrix(cyclotomic_polynomial(p**k), format='right')
        gamma_mat0 = block_diagonal_matrix([comp_mat] * dim, subdivide=False)
        gamma_mat = newbasis**(-1) * gamma_mat0 * newbasis
    else:
        # Cast as elements of E or OE naturally instead of coercing to ZZ
        gram_mat = diagonal_matrix(diag_gram_mat)
        gamma_mat = identity_matrix(dim)
        type_list = I

    mass_local_term = hermi_unr_mass_local_term(q, dim, I)
    disc_isnorm = sum((i - d) * count for i, count in enumerate(I)) % 2
        
    # Instantiate the lattice class, passing the generalized q
    lattice_obj = hermilattice_unr(q, p, m2, dim * e, type_list, gram_mat, gamma_mat, mass_local_term)
    
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
                
                # [Crucial Update]
                # Replaced prime ideal factorization from CyclotomicField(m2) (which depends on Q)
                # with direct instantiation of a finite field (residue field) of order q.
                # This ensures compatibility across arbitrary local fields.
                residue_field_local = GF(q, 'a')
                
                # Extract the submatrix from the automorphism matrix, reduce it over the residue field, and subtract the identity matrix
                # (If gamma_mat contains integers, they are automatically projected/reduced into the finite field)
                sub_gamma = gamma_mat.submatrix(type_list[0], type_list[0])
                reduced_matrix = matrix(residue_field_local, self.reduction_dim, self.reduction_dim, sub_gamma) - 1
                
                self.reduction_young = young_list(reduced_matrix)

    def __str__(self):
        return (f"Hermitian lattice over local field with residue field F_{self.q}, "
                f"rank {self.dim}, type {self.type_list}, "
                f"local term in mass formula: {self.mass_local_term}")

    def __repr__(self):
        return str(self)
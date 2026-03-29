load("gl_density.sage")
load("young.sage")

def build_gl_lattice(p, k, m2, d, e, dim, I, qv):
    """
    Construct a GL lattice object from a specified formal type.
    
    INPUT:
        p   : integer; a prime number
        k   : integer; exponent of p in the cyclotomic order
        m2  : integer; tame cyclotomic order
        d   : integer; degree parameter
        e   : integer; ramification index
        dim : integer; dimension/multiplicity
        I   : list of integers; type list representing the lattice invariants
        qv : integer; the size of the residue field of F_v (i.e., p^{fv})

    OUTPUT:
        A list `[disc_dummy, lattice_obj]`, where:
            - disc_dummy : a dummy invariant (Mod(0, 2))
            - lattice_obj: an instance of the gl_lattice class
    """
    if k > 0:
        comp_mat = companion_matrix(cyclotomic_polynomial(p**k), format='right')
        gamma_mat = block_diagonal_matrix([comp_mat] * dim, subdivide=False)
    else:
        gamma_mat = identity_matrix(dim)

    mass_local_term = gl_mass_local_term(qv, I)
    
    disc_dummy = Mod(0, 2) 
    
    lattice_obj = gl_lattice(p, m2, dim * e, gamma_mat, mass_local_term)
    return [disc_dummy, lattice_obj]


class gl_lattice:
    """
    Represents an unramified GL lattice over the local ring of integers.
    """
    def __init__(self, p, m2, dim, gamma_mat, mass_local_term):
        self.p = p
        self.m2 = m2
        self.dim = dim
        self.gamma_mat = gamma_mat  
        self.mass_local_term = mass_local_term

    def __repr__(self):
        return f"GL Lattice(mass={self.mass_local_term})"
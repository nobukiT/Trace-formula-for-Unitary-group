class QuadMod4:
    """
    Represents a unimodular quadratic lattice modulo 4.
    """
    def __init__(self, dim, evenpart, oddpart):
        """
        Initialize the QuadMod4 lattice.

        INPUT:
            dim      : integer; the total dimension of the lattice.
            evenpart : integer or Mod; the even part invariant (0 or 1 modulo 2).
            oddpart  : integer; the odd part invariant (0, 1, 2, 3, 4, or 5).
        """
        if dim < 0:
            raise ValueError("QuadMod4: negative dimension.")
        self.dim = dim
        self.evenpart = Mod(evenpart, 2)
        self.oddpart = oddpart
        
        # Determine the dimension of the odd part based on the oddpart invariant
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
        """
        Return the negation of the lattice modulo 4.
        """
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
    """
    Projects the basis onto the orthogonal complement of the plane spanned by indices i and j (modulo 4).
    This operation modifies the matrix in-place.

    INPUT:
        A : matrix; a symmetric matrix over Z/4Z.
        i : integer; first index of the plane basis.
        j : integer; second index of the plane basis.
        
    OUTPUT:
        None (the matrix A is modified in-place).
    """
    n = A.nrows()
    detinv = A[i, i] * A[j, j] - A[i, j] * A[j, i]
    for k in range(n):
        if k != i and k != j:
            a = detinv * (A[j, j] * A[k, i] - A[i, j] * A[k, j])
            b = detinv * (-A[j, i] * A[k, i] + A[i, i] * A[k, j])
            if a != 0 or b != 0:
                for l in range(n):
                    A[k, l] -= (a * A[i, l] + b * A[j, l])
                    A[l, k] = A[k, l] # Maintain symmetry
    return

def QuadMod4_from_mat(S):
    """
    Extracts the QuadMod4 invariant object from a given symmetric matrix.

    INPUT:
        S : matrix; a square symmetric matrix representing the lattice.
        
    OUTPUT:
        QuadMod4 object representing the unimodular invariants of S modulo 4.
    """
    dim = S.nrows()
    if dim != S.ncols():
        raise ValueError("Matrix must be square.")
    
    A = matrix(Integers(4), S)
    n = dim
    evenpart = Mod(0, 2)
    
    # Reduction loop
    while n > 2:
        # Adjustment to make A[0,0] even
        if Mod(A[0, 0], 2) != 0:
            if Mod(A[1, 1], 2) != 0:
                A.add_multiple_of_row(0, 1, 1)
                A.add_multiple_of_column(0, 1, 1)
            else:
                A.swap_rows(0, 1)
                A.swap_columns(0, 1)
        
        # Search for an even basis vector
        i = 1
        while i < n and (Mod(A[i, i], 2) != 0 or Mod(A[0, i], 2) == 0):
            i += 1
            
        if i == n: # If no even plane is found
            i = 1
            while i < n and (Mod(A[i, i], 2) == 0 or Mod(A[0, i], 2) == 0):
                i += 1
            if i == n: raise ValueError("Not unimodular.")
            
            j = 1
            while j < n and (Mod(A[j, j], 2) == 0 or Mod(A[0, j], 2) != 0):
                j += 1
            
            if j == n: # Orthogonal complement is even
                A.swap_rows(i, n - 1); A.swap_columns(i, n - 1)
                A.swap_rows(0, n - 2); A.swap_columns(0, n - 2)
                planeorthmod4_inplace(A, n - 2, n - 1)
                # Re-search
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

    # Final determination step
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
        
        # Does it represent 1?
        if A[0, 0] == 1 or A[1, 1] == 1 or (A[0, 0] + A[1, 1] + 2 * A[0, 1]) % 4 == 1:
            oddpart = 3 if A.det() == 1 else 4
        else:
            oddpart = 5
            
    return QuadMod4(dim, evenpart, oddpart)

def planeorthmod4_inplace_newbasis(A, n, P, i, j):
    """
    Orthogonalizes the rest of the basis with respect to the plane specified by indices i and j.
    Modifies both the matrix A and the basis transformation matrix P in-place.

    INPUT:
        A : matrix; a symmetric matrix over Z/4Z (modified in-place).
        n : integer; the dimension of the matrix.
        P : matrix; a basis transformation matrix over Z/2Z (modified in-place).
        i : integer; first index of the plane.
        j : integer; second index of the plane.
        
    OUTPUT:
        None
    """
    # Calculate the determinant (or its inverse) of the plane (i, j) over Z/4Z.
    # Assuming unimodularity, det is either 1 or 3 (=-1), so det == det_inv
    det_inv = A[i, i] * A[j, j] - A[i, j] * A[j, i]
    
    # Apply orthogonalization to all other basis vectors k
    for k in range(n):
        if k == i or k == j:
            continue
            
        # Find coefficients to project the k-th vector onto the plane.
        # a, b are such that A[k] = a*A[i] + b*A[j] (mod orthogonal complement of the plane)
        a = det_inv * (A[j, j] * A[k, i] - A[i, j] * A[k, j])
        b = det_inv * (-A[j, i] * A[k, i] + A[i, i] * A[k, j])
        
        # Subtract the i and j components from the k-th row/column to orthogonalize
        A.add_multiple_of_column(k, i, -a)
        A.add_multiple_of_column(k, j, -b)
        A.add_multiple_of_row(k, i, -a)
        A.add_multiple_of_row(k, j, -b)
        
        # Reflect this operation on the transformation matrix P (over Z/2Z)
        if Mod(a, 2) != 0:
            P.add_multiple_of_column(k, i, 1)
        if Mod(b, 2) != 0:
            P.add_multiple_of_column(k, j, 1)

def QuadMod4_newbasis_from_mat(S):
    """
    Converts a unimodular lattice S into normal form over Z/4Z and returns 
    its invariants along with the basis transformation matrix.

    INPUT:
        S : matrix; a square symmetric matrix.

    OUTPUT:
        tuple (QuadMod4, Matrix); the invariants object and the transformation matrix over Z/2Z.
    """
    dim = S.nrows()
    if dim != S.ncols():
        raise ValueError("Error in quadmod4_newbasis: nonsquare matrix cannot be symmetric.")

    # Set up rings for computation
    R4 = Integers(4)
    R2 = Integers(2)
    
    A = matrix(R4, S)
    new_basis = identity_matrix(R2, dim)
    n = dim
    even_part = Mod(0, 2)

    # 1. Loop to extract even unimodular planes
    while n > 2 or (n == 2 and Mod(A[0, 0], 2) == 0 and Mod(A[1, 1], 2) == 0):
        # Make the last basis vector A[n-1] even
        if Mod(A[n-1, n-1], 2) != 0:
            if Mod(A[n-2, n-2], 2) != 0:
                # Replace e_n with e_n + e_{n-1}
                A[n-1, n-1] += 2 * A[n-1, n-2] + A[n-2, n-2]
                for i in range(n - 1):
                    val = A[n-2, i]
                    A[n-1, i] += val
                    A[i, n-1] += val
                new_basis.add_multiple_of_column(n-1, n-2, 1)
            else:
                # Swap e_{n-1} and e_n
                A.swap_columns(n-2, n-1)
                A.swap_rows(n-2, n-1)
                new_basis.swap_columns(n-2, n-1)

        # Search for an even basis vector orthogonal to e_n
        idx = -1
        for i in range(n - 2, -1, -1):
            if Mod(A[i, i], 2) == 0 and Mod(A[n-1, i], 2) != 0:
                idx = i
                break
        
        if idx < 0:
            # If no even basis vector is found, search for an odd one
            for i in range(n - 2, -1, -1):
                if Mod(A[i, i], 2) != 0 and Mod(A[n-1, i], 2) != 0:
                    idx = i
                    break
            
            if idx < 0:
                raise RuntimeError("Error: quadmod4_newbasis: not unimodular.")

            # Search for an odd basis vector orthogonal to e_n
            j_idx = -1
            for j in range(n - 2, -1, -1):
                if Mod(A[j, j], 2) != 0 and Mod(A[n-1, j], 2) == 0:
                    j_idx = j
                    break
            
            if j_idx < 0:
                # (e_i, e_n) forms an odd unimodular plane. Make the orthogonal complement even.
                if idx > 0:
                    A.swap_columns(idx, 0); A.swap_rows(idx, 0)
                    new_basis.swap_columns(idx, 0)
                A.swap_columns(1, n-1); A.swap_rows(1, n-1)
                new_basis.swap_columns(1, n-1)
                
                # Plane orthogonalization
                planeorthmod4_inplace_newbasis(A, n, new_basis, 0, 1)
                
                # Re-search
                for i in range(n - 2, -1, -1):
                    if Mod(A[i, i], 2) == 0 and Mod(A[n-1, i], 2) != 0:
                        idx = i
                        break
                if idx < 0:
                    raise RuntimeError("Error: quadmod4_newbasis reduction failed.")
            else:
                # (e_i + e_j, e_n) will form an even unimodular plane
                A.add_multiple_of_column(idx, j_idx, 1)
                A.add_multiple_of_row(idx, j_idx, 1)
                new_basis.add_multiple_of_column(idx, j_idx, 1)

        # Finalize (e_idx, e_{n-1}) as an even unimodular plane
        if idx < n - 2:
            A.swap_columns(n-2, idx); A.swap_rows(n-2, idx)
            new_basis.swap_columns(n-2, idx)
        
        planeorthmod4_inplace_newbasis(A, n, new_basis, n-2, n-1)

        # Adjustment for standard forms diag(0,2) or diag(2,0)
        if A[n-2, n-2] == 0 and A[n-1, n-1] == 2:
            A.add_multiple_of_column(n-1, n-2, 1); A.add_multiple_of_row(n-1, n-2, 1)
            new_basis.add_multiple_of_column(n-1, n-2, 1)
        elif A[n-2, n-2] == 2 and A[n-1, n-1] == 0:
            A.add_multiple_of_column(n-2, n-1, 1); A.add_multiple_of_row(n-2, n-1, 1)
            new_basis.add_multiple_of_column(n-2, n-1, 1)

        # Update the even_part invariant
        if even_part == 0:
            if A[n-1, n-1] == 2:
                even_part = Mod(1, 2)
        else:
            if A[n-1, n-1] == 0:
                # Swap planes
                indices = [n-2, n-1, n, n+1]
                for k in range(0, 2):
                    A.swap_rows(indices[k], indices[k+2])
                    A.swap_columns(indices[k], indices[k+2])
                    new_basis.swap_columns(indices[k], indices[k+2])
            else:
                even_part = Mod(0, 2)
                # Clear bases and update
                for k in [n-2, n-1, n, n+1]: A[k, k] = 0
                new_basis.add_multiple_of_column(n-2, n, 1)
                new_basis.add_multiple_of_column(n+1, n-1, 1)
                new_basis.add_multiple_of_column(n-1, n-2, 1)
                new_basis.add_multiple_of_column(n, n+1, 1)
        
        n -= 2

    # 2. Classification of the remaining odd part
    odd_part = 0
    if n == 1:
        odd_part = 1 if A[0, 0] == 1 else 2
    elif n == 2:
        # Adjustment when diagonal entries are 0
        if Mod(A[0, 0], 2) == 0:
            A.add_multiple_of_column(0, 1, 1); A.add_multiple_of_row(0, 1, 1)
            new_basis.add_multiple_of_column(0, 1, 1)
        elif Mod(A[1, 1], 2) == 0:
            A.add_multiple_of_column(1, 0, 1); A.add_multiple_of_row(1, 0, 1)
            new_basis.add_multiple_of_column(1, 0, 1)

        # Consolidation of special cases
        if even_part == 1 and A[0, 0] * A[1, 1] == 1:
            even_part = 0
            # Operations to resolve interactions between A[0,0], A[1,1] and A[2,2], A[3,3]
            for i in [0, 1]:
                for j in [2, 3]:
                    A.add_multiple_of_column(i, j, 1)
                    A.add_multiple_of_column(j, i, 1)
                    A.add_multiple_of_row(i, j, 1)
                    A.add_multiple_of_row(j, i, 1)
                    new_basis.add_multiple_of_column(i, j, 1)
                    new_basis.add_multiple_of_column(j, i, 1)

        # Determine the invariant of the odd part
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
    """
    Represents a regular quadratic form over the p-adic field Q_p.
    """
    def __init__(self, p, dim, disc, hasse):
        """
        Initialize the QuadQp structure.

        INPUT:
            p     : integer; a prime number.
            dim   : integer; the dimension of the quadratic form.
            disc  : integer or rational; the discriminant.
            hasse : integer (+1 or -1); the Hasse invariant.
        """
        if not is_prime(p): raise ValueError(f"{p} is not prime")
        if dim <= 0: raise ValueError("Dimension must be >= 1")
        if disc == 0 or hasse**2 != 1: raise ValueError("Invalid disc or Hasse invariant")
        
        self.p = p
        self.dim = dim
        self.disc = self._unique_repr_mod_squares(disc, p)
        self.hasse = hasse
        
        # Dimension consistency checks
        if dim == 1 and hasse == -1:
            raise ValueError("Impossible form in dimension 1")
        if dim == 2 and hasse == -1:
            if (p > 2 and kronecker(-self.disc, p) == 1) or (p == 2 and Mod(self.disc, 8) == 7):
                raise ValueError("Impossible form in dimension 2")

    @staticmethod
    def _unique_repr_mod_squares(d, p):
        """
        Returns a unique representative for the discriminant modulo squares in Q_p^x.
        """
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
        """
        Checks whether the discriminant is equivalent to d modulo squares.
        """
        v = d.valuation(self.p)
        if (v - self.disc.valuation(self.p)) % 2 != 0:
            return False
        d_unit = d // self.p**v if v % 2 == 0 else d // self.p**(v - 1)
        if self.p > 2:
            return kronecker(self.disc / d_unit, self.p) == 1
        else:
            return Mod(self.disc / d_unit, 8) == 1

    def __add__(self, other):
        """
        Returns the orthogonal direct sum of two p-adic quadratic forms.
        """
        if self.p != other.p: raise ValueError("Primes do not match.")
        new_hasse = self.hasse * other.hasse * hilbert_symbol(self.disc, other.disc, self.p)
        return QuadQp(self.p, self.dim + other.dim, self.disc * other.disc, new_hasse)

    def contains_regular_lattice(self):
        """
        Checks if the quadratic space contains a regular integral lattice.
        
        OUTPUT:
            Boolean.
        """
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
    """
    Constructs a QuadQp object safely. A is guaranteed to be over QQ.
    """
    dim = A.nrows()
    A_QQ = matrix(QQ, A)
    disc = A_QQ.det()
    
    if disc == 0:
        raise ValueError(f"Discriminant is zero at prime {p}. Matrix:\n{A_QQ}")
    
    Q = QuadraticForm(QQ, 2 * A_QQ)
    hasse = Q.hasse_invariant(p)
    
    return QuadQp(p, dim, disc, 1 if hasse > 0 else -1)

def sum_list_quad(L):
    """
    Computes the orthogonal direct sum of a list of quadratic forms.

    INPUT:
        L : list; a list of QuadQp objects.
        
    OUTPUT:
        QuadQp object representing the total sum.
    """
    if not L: raise ValueError("Empty list.")
    res = L[0]
    for item in L[1:]:
        res += item
    return res

class QuadGlob:
    """
    Represents a non-degenerate quadratic form over the rational numbers QQ.
    """
    def __init__(self, dim, disc, conductor, negdim):
        """
        Initialize the QuadGlob form.

        INPUT:
            dim       : integer; the total dimension of the quadratic form.
            disc      : integer; the discriminant (will be reduced to square-free).
            conductor : integer; the conductor.
            negdim    : integer; the dimension of the negative-definite subspace (signature info).
        """
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
        """
        Reduces an integer modulo squares (returns the square-free part).
        """
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
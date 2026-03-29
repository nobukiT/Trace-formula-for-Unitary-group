def get_ref_bil_form_munip_3(m):
    """
    Computes the reference bilinear form matrix S for a Jordan block of size 2m 
    associated with a minus unipotent isometry over F_3.

    INPUT:
        m : integer; half the size of the target Jordan block (size = 2m).
        
    OUTPUT:
        matrix: A 2m x 2m symmetric matrix over Integers(3) representing the reference bilinear form.
    """
    if m < 1:
        return matrix(Integers(3), 0, 0)
        
    # Prepare the polynomial ring and the companion matrix
    r0 = PolynomialRing(Integers(3), 'x0')
    x0 = r0.gen()
    poly = (x0 + 1)**(2 * m)
    p_mat = companion_matrix(poly)
    
    # Computation in the quotient ring
    r_quot = r0.quotient(poly, 'x')
    x = r_quot.gen()
    
    dim = 2 * m
    s_mat = matrix(Integers(3), dim, dim)
    for i in range(dim):
        for j in range(dim):
            # Get the list of coefficients for x^(j-i+m)
            coeffs = (x**(j - i + m)).list()
            # Safely extract the required coefficients
            c0 = coeffs[0] if len(coeffs) > 0 else 0
            cn = coeffs[dim - 1] if len(coeffs) >= dim else 0
            s_mat[i, j] = c0 + cn
            
    # Construct a new basis
    e1 = matrix(Integers(3), dim, 1, [1] + [0] * (dim - 1))
    basis_cols = [e1]
    for i in range(dim - 1):
        # Apply P to the previous column vector to generate the next one
        basis_cols.append(basis_cols[i] + p_mat * basis_cols[i])
        
    new_basis = block_matrix(1, dim, basis_cols)
    # Return new_s = transpose(basis) * S * basis
    return new_basis.transpose() * s_mat * new_basis


def invariants_munip_bil_form_3(s_mat):
    """
    Calculates the invariants of a 'minus unipotent' bilinear form over F_3.

    INPUT:
        s_mat : matrix; a square matrix over Integers(3) representing the bilinear form.
        
    OUTPUT:
        tuple: A sequence of invariants where the i-th entry corresponds to Jordan blocks 
               of size (i+1). Odd sizes contain just the multiplicity, while even sizes 
               contain a tuple (multiplicity, discriminant mod 2).
    """
    if s_mat.nrows() == 0:
        return ()
        
    # Calculate P = transpose(S^-1) * S
    p_mat = (s_mat.inverse().transpose()) * s_mat
    n_mat = 1 + p_mat  # Nilpotent part
    
    # Compute the Jordan normal form and the transformation matrix
    new_n, new_basis = n_mat.jordan_form(transformation=True)
    # Extract block sizes and their multiplicities
    young_tableau = get_jordan_mults(new_n.subdivisions()[0], n_mat.nrows())
    
    # Normalize the basis (adjustment using an anti-diagonal matrix)
    block_list = [antidiag_scalar_matrix(Integers(3), r[0], 1) for r in young_tableau for _ in range(r[1])]
    new_basis *= block_diagonal_matrix(block_list, subdivide=False)
    
    # The bilinear form S after transformation
    transformed_s = new_basis.transpose() * s_mat * new_basis
    
    res_list = [0] * young_tableau[0][0]
    cur_pos = 0
    
    for r_dim, r_mult in young_tableau:
        if r_dim % 2 == 0:
            # For even-sized Jordan cells, extract the corresponding submatrix to determine the discriminant.
            # Call the reference matrix generator get_ref_bil_form_munip_3(m) on the fly.
            ref_s_inv = get_ref_bil_form_munip_3(r_dim // 2).inverse()
            q_block_inv = block_diagonal_matrix([ref_s_inv] * r_mult, subdivide=False)
            
            # Extract the submatrix
            s_sub = transformed_s[cur_pos:cur_pos + r_dim * r_mult, cur_pos:cur_pos + r_dim * r_mult]
            q_sub = q_block_inv * s_sub
            
            # Calculate the reduced discriminant by extracting representative values
            q_red = matrix(Integers(3), r_mult, r_mult, lambda a, b: q_sub[a * r_dim, b * r_dim])
            disc = Mod(0 if kronecker(q_red.det(), 3) == 1 else 1, 2)
            res_list[r_dim - 1] = (r_mult, disc)
        else:
            # For odd sizes, only the multiplicity is recorded
            res_list[r_dim - 1] = r_mult
        cur_pos += r_dim * r_mult
        
    # Format the result (fill in empty slots for even sizes with zero multiplicity)
    for i in range(len(res_list)):
        if (i + 1) % 2 == 0 and res_list[i] == 0:
            res_list[i] = (0, Mod(0, 2))
            
    return tuple(res_list)


def antidiag_scalar_matrix(base_ring, dim, scalar):
    """
    Generates a matrix filled with a specified scalar along its anti-diagonal.

    INPUT:
        base_ring : Ring; the base ring of the matrix (e.g., Integers(3)).
        dim       : integer; the dimension of the square matrix.
        scalar    : element; the value to place on the anti-diagonal.
        
    OUTPUT:
        matrix: A dim x dim anti-diagonal matrix.
    """
    a = matrix(base_ring, dim, dim)
    for i in range(dim):
        a[i, dim - i - 1] = scalar
    return a


def get_jordan_mults(subdivisions, total_dim):
    """
    Calculates the block sizes and their multiplicities from Jordan decomposition subdivisions.

    INPUT:
        subdivisions : list; a list of partition indices from the Jordan form subdivisions.
        total_dim    : integer; the total dimension of the matrix.
        
    OUTPUT:
        list: A list of pairs [block_size, multiplicity].
    """
    if not subdivisions:
        return [[total_dim, 1]]
        
    res = []
    prev_pos = 0
    current_dim = 0
    current_mult = 0
    
    for i in subdivisions:
        dim_diff = i - prev_pos
        if dim_diff == current_dim:
            current_mult += 1
        else:
            if current_dim > 0:
                res.append([current_dim, current_mult])
            current_dim = dim_diff
            current_mult = 1
        prev_pos = i
        
    # Process the final block
    last_dim = total_dim - prev_pos
    if last_dim == current_dim:
        res.append([current_dim, current_mult + 1])
    else:
        if current_mult > 0:
            res.append([current_dim, current_mult])
        res.append([last_dim, 1])
    return res


def core_quad_3(gram_mat):
    """
    Calculates the dimension and the quadratic residuosity of the core (non-degenerate part) 
    of a quadratic form over F_3.

    INPUT:
        gram_mat : matrix; a symmetric matrix over Integers(3).
        
    OUTPUT:
        tuple: (rank, discriminant modulo 2), where discriminant is 0 if it is a square, and 1 otherwise.
    """
    # Obtain the indices of the non-degenerate subspace
    kern = gram_mat.right_kernel().echelonized_basis()
    # In practice, representative indices of the orthogonal complement are extracted from the echelonized basis
    # (Maintaining existing logic)
    dim = gram_mat.nrows()
    rank = gram_mat.rank()
    
    # Simplified discriminant determination
    # (While a full inverse might theoretically be needed, this logic computes the det of the reduced matrix)
    if rank == 0:
        return (0, Mod(0, 2))
        
    # Extract non-zero indices from the echelon form
    pivot_indices = gram_mat.pivots()
    sub_det = gram_mat[pivot_indices, pivot_indices].det()
    
    disc = Mod(0 if kronecker(sub_det, 3) == 1 else 1, 2)
    return (len(pivot_indices), disc)


def invariants_quad_with_unip_3(gram_mat, gamma_mat):
    """
    Calculates the invariants of a quadratic form over F_3 equipped with a unipotent action.

    INPUT:
        gram_mat  : matrix; the Gram matrix of the quadratic form over Integers(3).
        gamma_mat : matrix; the automorphism (isometry) matrix over Integers(3).
        
    OUTPUT:
        tuple: (total_dimension, total_discriminant mod 2, detailed_invariants)
               where detailed_invariants is a tuple starting with the core invariants on the 
               kernel of (1 - gamma), followed by the invariants of the minus unipotent bilinear form.
    """
    t_mat = 1 - gamma_mat
    dim_tot = gram_mat.nrows()
    
    # Calculate the subspace Ker(1 - gamma)
    kern_basis = t_mat.right_kernel().basis_matrix()
    s_on_ker = kern_basis * gram_mat * kern_basis.transpose()
    
    # Extract the part associated with the non-trivial action
    pivots = t_mat.pivots()
    s_mat = (gram_mat * t_mat)[pivots, pivots]
    
    disc_sq = Mod(0 if kronecker(gram_mat.det(), 3) == 1 else 1, 2)
    inv_bil = invariants_munip_bil_form_3(s_mat)
    
    return (dim_tot, disc_sq, (core_quad_3(s_on_ker),) + inv_bil)
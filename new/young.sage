def young_list(M):
    """
    Calculate the size of Jordan blocks (Young diagram) for a nilpotent matrix M.
    
    INPUT:
        M -- a square matrix (must be nilpotent)

    OUTPUT:
        List of integers representing the sizes of the Jordan blocks (a partition).
    """
    if not M.is_square():
        raise ValueError("The matrix must be a square matrix.")
        
    total_dim = M.nrows()
    if total_dim == 0:
        return []

    partition = []
    current_nullity = total_dim - M.rank()  # dim(ker(M))
    prev_nullity = 0
    
    # Variable to keep track of the matrix powers
    M_pow = M
    
    # Loop until the dimension of the kernel stops increasing
    while current_nullity > prev_nullity:
        # This increment forms the length of a "column" in the Young diagram
        partition.append(current_nullity - prev_nullity)
        
        # Compute the nullity for the next power
        prev_nullity = current_nullity
        M_pow *= M
        current_nullity = total_dim - M_pow.rank()

    return partition
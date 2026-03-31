def young_list(M):
    """
    Return the sizes of the Jordan blocks of a nilpotent matrix M.

    INPUT:
        M -- a square nilpotent matrix

    OUTPUT:
        List of integers representing the sizes of the Jordan blocks,
        in decreasing order.
    """
    if not M.is_square():
        raise ValueError("The matrix must be a square matrix.")
        
    n = M.nrows()
    if n == 0:
        return []

    columns = []
    prev_nullity = 0
    M_pow = M

    while True:
        current_nullity = n - M_pow.rank()
        if current_nullity == prev_nullity:
            break
        columns.append(current_nullity - prev_nullity)
        prev_nullity = current_nullity
        M_pow *= M

    blocks = []
    k = 1
    while True:
        count = sum(1 for c in columns if c >= k)
        if count == 0:
            break
        blocks.append(count)
        k += 1

    return blocks

load("utils.sage")

def relative_L_special_value(Ei, Fi, s, prec=1000, verbose=False):
    """
    Return L(E/F, 1-j) for odd j >= 1.
    Evaluates the exact limit avoiding 0/0 by computing the d-th derivatives using PARI.
    """
    j = 1 - s

    if j % 2 == 0 or j < 1:
        raise ValueError(f"j must be odd positive, got j={j}")

    if j == 1:
        return QQ(-Ei.class_number()) / QQ(Ei.number_of_roots_of_unity())

    d = Fi.degree()

    R = RealField(prec)
    s_val = R(s)

    if verbose:
        print(f"Evaluating exact limit at s={s} using {d}-th derivatives.")
        print(f"   prec = {prec} bits")

    nf_E = Ei.pari_nf()
    
    if Fi == QQ:
        val_E_deriv = pari.lfun(nf_E, s_val, d)
        val_F_deriv = pari.lfun(1, s_val, d) 
    else:
        nf_F = Fi.pari_nf()
        val_E_deriv = pari.lfun(nf_E, s_val, d)
        val_F_deriv = pari.lfun(nf_F, s_val, d)

    if val_F_deriv == 0:
        raise ZeroDivisionError(f"d-th derivative of zeta_F is zero at s={s}")

    approx_sage = R((val_E_deriv / val_F_deriv).real())
    q = exactify_to_QQ(approx_sage)

    if verbose:
        print(f"   approx L(E/F,{s}) = {approx_sage}")
        print(f"   recovered         = {q}")

    return q


def unitary_block_L_product(block_size, Ei, Fi, prec=200, verbose=False):
    """
    Compute

        prod_{j=1}^{block_size} term_j

    where
        term_j = L(E/F, 1-j)   if j is odd,
               = zeta_F(1-j)   if j is even.

    Returns an exact rational whenever possible.
    """
    ans = QQ(1)
    zF = Fi.zeta_function()

    for j in range(1, block_size + 1):
        s = 1 - j

        if j % 2 == 0:
            # even j -> trivial character
            raw_term = zF(s)
            term = exactify_to_QQ(raw_term)

            if verbose:
                print(f"j={j}, s={s}: zeta_F({s}) = {term}")

        else:
            # odd j -> relative quadratic character
            term = relative_L_special_value(Ei, Fi, s, prec=prec, verbose=verbose)

            if verbose:
                print(f"j={j}, s={s}: L(E/F, {s}) = {term}")

        ans *= term

    return QQ(ans)

load("utils.sage")

def relative_L_special_value(Ei, Fi, s, prec=1000, verbose=False):
    """
    Return L(E/F, 1-j) for odd j >= 1, assuming E/F is a CM extension.

    Here
        L(E/F, s) = zeta_E(s) / zeta_F(s).

    For j > 1, this evaluates the exact limit at s = 1-j by taking
    the d-th derivatives (where d = [F:Q]) of zeta_E and zeta_F.

    For j = 1 (i.e. s = 0), it uses the exact class number formula:
        L(E/F, 0) = (2^{d-1} / Q) * (h_E / h_F) * (w_F / w_E),
    where d = [F:Q], h_K is the class number, w_K is the number
    of roots of unity, and Q is the Hasse unit index.
    """
    j = 1 - s

    if j % 2 == 0 or j < 1:
        raise ValueError(f"j must be odd positive, got j={j}")

    d = Fi.degree() if Fi != QQ else 1

    if j == 1:
        hE = QQ(Ei.class_number())
        wE = QQ(Ei.number_of_roots_of_unity())
        RE = Ei.regulator()

        if Fi == QQ:
            hF = QQ(1)
            wF = QQ(2)
            RF = QQ(1)
        else:
            hF = QQ(Fi.class_number())
            wF = QQ(Fi.number_of_roots_of_unity())
            RF = Fi.regulator()
       

        if verbose:
            print(f"Exact class number formula at s=0:")
            print(f"  d={d}, hE={hE}, hF={hF}, wE={wE}, wF={wF}, Q={Q}")

        exact_val = (RE/RF) * (hE / hF) * (wF / wE)
        return exact_val

    R = RealField(prec)
    s_val = R(s)

    if verbose:
        print(f"Evaluating exact limit at s={s} using {d}-th derivatives.")
        print(f"  prec = {prec} bits")

    nf_E = Ei.pari_nf()
    deriv_order = d 

    if Fi == QQ:
        val_E_deriv = pari.lfun(nf_E, s_val, deriv_order)
        val_F_deriv = pari.lfun(1, s_val, deriv_order)
    else:
        nf_F = Fi.pari_nf()
        val_E_deriv = pari.lfun(nf_E, s_val, deriv_order)
        val_F_deriv = pari.lfun(nf_F, s_val, deriv_order)

    if val_F_deriv == 0:
        raise ZeroDivisionError(f"{deriv_order}-th derivative of zeta_F is zero at s={s}")

    val = val_E_deriv / val_F_deriv
    if abs(val.imag()) > 2**(-prec//3):
        raise ValueError(f"Unexpected imaginary part: {val.imag()}")

    approx_sage = R(val.real())
    q = exactify_to_QQ(approx_sage)

    if verbose:
        print(f"  approx L(E/F,{s}) = {approx_sage}")
        print(f"  recovered         = {q}")

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
    pari_Fi = pari(Fi)

    for j in range(1, block_size + 1):
        s = 1 - j

        if j % 2 == 0:
            # even j -> trivial character
            raw_term = pari.lfun(pari_Fi, s)
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
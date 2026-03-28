

def block_zeta_L_data_unitary(m, D_K, verbose=False):
    """
    Return the block fields Ei/Fi together with zeta and Hecke-L data.
    OUTPUT:
        dict with keys:
            "m"
            "Ei"
            "Fi"
            "zeta_Fi"       -- Dedekind zeta function of Fi
            "zeta_Ei"       -- Dedekind zeta function of Ei
            "L_Ei_over_Fi"  -- represented as zeta_Ei / zeta_Fi
    """
    data = block_fields_unitary(m, D_K, verbose=verbose)
    Ei = data["Ei"]
    Fi = data["Fi"]
    data["zeta_Ei"] = Ei.zeta_function()
    data["zeta_Fi"] = Fi.zeta_function()
    data["L_Ei_over_Fi"] = data["zeta_Ei"] / data["zeta_Fi"]
    return data

def unitary_block_L_product(block_size, Ei, Fi, verbose=False):
    """
    Return the exact special value

        prod_{j=1}^{block_size} L_Fi(1-j, chi_{Ei/Fi}^j)

    OUTPUT:
        an exact Sage expression whenever possible
    """
    ans = QQ(1)

    zeta_E = Ei.zeta_function()
    zeta_F = Fi.zeta_function()

    for j in range(1, block_size + 1):
        s = 1 - j

        if j % 2 == 0:
            # trivial character
            term = zeta_F(s)
            if verbose:
                print(f"j={j}, s={s}: zeta_F({s}) = {term}")
        else:
            # quadratic character
            term = zeta_E(s) / zeta_F(s)
            if verbose:
                print(f"j={j}, s={s}: zeta_E({s}) / zeta_F({s}) = {term}")
        ans *= term
    return ans
    

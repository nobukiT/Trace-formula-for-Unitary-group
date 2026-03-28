load("utils.sage")

def block_fields_unitary(m, D_K, verbose=False):
    """
    construct the corresponding fields E_i/F_i for the centralizer block.
    OUTPUT:
        dict with keys:
            "m"          : m
            "E_base"     : the quadratic field E = Q(sqrt(D_K))
            "cyclo"      : Q(zeta_m)
            "Ei"         : compositum E(zeta_m)
            "Fi"         : maximal totally real subfield (when Ei is CM)
            "split"      : whether Ei is not a field-CM block in the expected sense
            "deg_Ei"     : [Ei:Q]
            "deg_Fi"     : [Fi:Q]
    """
    QQx.<x> = QQ[]
    E0.<w> = QuadraticField(D_K)

    # m=1,2 are degenerate but still valid
    if m == 1:
        L = QQ
    else:
        L = CyclotomicField(m)

    # compositum E(zeta_m)
    Ei, emb1, emb2 = E0.composite_fields(L, both_maps=True)[0]

    out = {
        "m": m,
        "E_base": E0,
        "cyclo": L,
        "Ei": Ei,
        "deg_Ei": Ei.degree(),
    }
    Fi = Ei.maximal_totally_real_subfield()[0]
    out["Fi"] = Fi
    out["deg_Fi"] = Fi.degree()
    out["split"] = False

    return out
#============================================================================
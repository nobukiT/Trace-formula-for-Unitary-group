load("global_herm.sage")
load("orb_int_herm.sage")

def run_mass_test():
    orb_int_unitary_db.clear()
    hermi_latt_db.clear()
    gl_latt_db.clear()
    # --- テストパラメータの設定 ---
    n = 2
    N = 2 * n
    D_K = -3
    
    primes_list = [2, 3, 5, 7, 11, 13]  

    print("=========================================================")
    print(f" Step 1: 軌道積分（局所密度）データベースの構築")
    print(f" 群: U({n},{n}), 判別式 D_K: {D_K}")
    print("=========================================================")
    update_orb_int_unitary_db(n, primes_list, D_K, pretend=False)
    
    print("\n=========================================================")
    print(f" Step 2: 大域的な質量 (Mass Formula) の計算")
    print("=========================================================")
    results = mass_list_unitary(n, primes_list, orb_int_unitary_db, D_K)
    
    # --- 結果の集計 ---
    print("\n=========================================================")
    print(" 最終計算結果")
    print("=========================================================")
    total_mass = 0
    for conj_class, mass in results:
        # conj_class は [m, 重複度] のリストの集まり
        class_str = " * ".join([f"Phi_{m}(x)^{mult}" for m, mult in conj_class])
        print(f"Mass of [ {class_str} ]: {mass}")
        total_mass += mass

# テストの実行
if __name__ == "__main__":
    run_mass_test()
# ==========================================================
# UNITARY GROUP U(n) MASS FORMULA - INTEGRATION TEST
# ==========================================================
import time

load("utils.sage")
load("Lfunc.sage")
load("jordan.sage")
load("quadclass.sage")
load("quadlatticeclass.sage")
load("hermilatticeclass.sage")
load("hermi.sage")
load("gl.sage")
load("orb_int_herm.sage")
load("global_herm.sage")

def run_final_integration_test():
    """
    U(2) における類数（質量）計算の統合テストを実行します。
    """
    print("=" * 80)
    print("🚀 STARTING FINAL INTEGRATION TEST")
    print("=" * 80)
    
    n = 4
    D_K = -4  # E = Q(sqrt(-3))
    
    primes_list = list(primes(n+2))
    start_time = time.time()

    try:
        print("\n>>> TEST 1: update_orb_int_unitary_db")
        update_orb_int_unitary_db(n, primes_list, D_K, pretend=False, verbose=True)
        
        print("\nTEST 1 PASSED: Database updated successfully.")
        print("\n>>> TEST 2: mass_list_unitary")
        # global_herm.sage 内の関数を呼び出し
        # imposed_negdim=-1 はデフォルトの U(n/2, n/2) を計算
        results = mass_list_unitary(
            n, 
            primes_list, 
            orb_int_unitary_db, 
            D_K, 
            imposed_negdim=-1, 
            prec=100, 
            verbose=True
        )
        
        print("\nTEST 2 PASSED: Mass computation completed successfully.")
        
        # --------------------------------------------------
        # 結果の出力
        # --------------------------------------------------
        print("-" * 80)
        print("RESULTS (Stable Conjugacy Classes and their Masses):")
        if not results:
            print("  No classes found with non-zero mass.")
        for conj_class, mass in results:
            print(f"  Class: {conj_class}  ==>  Mass: {mass}")
        
        elapsed = time.time() - start_time
        print("-" * 80)
        print(f"Test finished in {elapsed:.2f} seconds.")
        print("=" * 80)

    except Exception as e:
        print(f"\n❌ INTEGRATION TEST FAILED")
        import traceback
        traceback.print_exc()
        print("=" * 80)

# テストの実行
if __name__ == "__main__":
    run_final_integration_test()
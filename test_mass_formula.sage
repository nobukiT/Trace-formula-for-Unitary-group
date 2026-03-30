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
    print("=" * 80)
    print("🚀 STARTING FINAL INTEGRATION TEST (Class Number 1 Fields)")
    print("=" * 80)
    
    n_list = [4]
    class_number_one_discriminants = [-3, -4, -7, -8]
    
    for n in n_list:
        print(f"\n\n" + "#" * 60)
        print(f"### TESTING FOR n = {n} ###")
        print("#" * 60)
        for D_K in class_number_one_discriminants:
            
            base_primes = list(primes(n + 2))
            ramified_primes = prime_divisors(abs(D_K))
            primes_list = sorted(list(set(base_primes + ramified_primes)))

            print(f"\n\n" + "#" * 60)
            print(f"### TESTING FOR D_K = {D_K} ###")
            print("#" * 60)

            try:
                print(f"\n>>> [D_K={D_K}] TEST 1: update_orb_int_unitary_db")
                update_orb_int_unitary_db(n, primes_list, D_K, pretend=False, verbose=True)

                print(f"\nTEST 1 PASSED: Database updated for D_K={D_K}.")
                print(f"\n>>> [D_K={D_K}] TEST 2: mass_list_unitary")

                # mass_list_unitary の呼び出し
                results = mass_list_unitary(
                    n, 
                    primes_list, 
                    orb_int_unitary_db, 
                    D_K, 
                    imposed_negdim=-1, 
                    prec=100, 
                    verbose=True
                )

                print(f"\nTEST 2 PASSED: Mass computation completed for D_K={D_K}.")

                # --------------------------------------------------
                print("-" * 60)
                print(f"RESULTS for D_K = {D_K}:")
                if not results:
                    print("  No classes found with non-zero mass.")
                else:
                    for conj_class, mass in results:
                        print(f"  Class: {conj_class}  ==>  Mass: {mass}")
                print("-" * 60)

            except Exception as e:
                print(f"\n❌ INTEGRATION TEST FAILED FOR D_K = {D_K}")
                import traceback
                traceback.print_exc()

        print("\n" + "=" * 80)
        print("✅ ALL TESTS FINISHED")
        print("=" * 80)

if __name__ == "__main__":
    run_final_integration_test()
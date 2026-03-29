# Unitary Group Mass Formula and Lattice Classification

A SageMath-based library for computing the **mass formula** (class numbers) of stable conjugacy classes of unitary groups \( U(p,q) \) over algebraic number fields.

This project provides an end-to-end computational pipeline, including:

- **local lattice classification**
- **local density and orbital integral computation**
- **global mass formula assembly**
- **exact special value evaluation of relative \( L \)-functions**

It is designed for arithmetic applications involving Hermitian lattices, quadratic forms, and automorphic mass computations.

---

## Mathematical Background

This implementation is inspired in part by the computational framework developed for split classical groups in the following work:

> Olivier Taïbi,  
> **Dimensions of spaces of level one automorphic forms for split classical groups using the trace formula**,  
> *Annales scientifiques de l'École Normale Supérieure* **50** (2017), no. 2, 269–344.  
> DOI: [10.24033/asens.2323](https://doi.org/10.24033/asens.2323)

In particular, the present codebase adapts and extends related ideas in the setting of **unitary groups**, with emphasis on:

- stable conjugacy class enumeration,
- local orbital integral calculations,
- lattice classification over local fields,
- and global mass formula computations.

It also incorporates additional routines specific to Hermitian lattices and ramified local structures that do not appear in the split orthogonal/symplectic setting.

If relevant, one may also compare with the tables and source code accompanying Taïbi’s paper.

---

## Features

### Mass Formula Computation
- Computes the total mass for stable conjugacy classes of unitary groups.
- Incorporates **Tamagawa factors** and **global \( L \)-function contributions**.

### Local Density and Orbital Integrals
- Computes local densities for Hermitian lattices in:
  - **unramified cases**
  - **ramified odd residue characteristic cases**
  - **the delicate \( p = 2 \) case**
- Supports **GL-type local density computations** for split primes.
- Includes a **database system** for precomputing and caching local orbital integral data.

### Lattice Classification
- Enumerates isomorphism classes of **Hermitian lattices over local fields**.
- Provides advanced tools for **quadratic lattices over \( \mathbb{Z}_p \)**, including specialized handling for:
  - \( p = 2 \)
  - \( p = 3 \)

### Invariant Theory and Local Structure
- Implements **generalized Jordan decomposition** for Hermitian and symmetric forms.
- Classifies **unimodular quadratic forms modulo 4**.
- Computes refined invariants for **quadratic forms with unipotent actions over \( \mathbb{F}_3 \)**.

### Number-Theoretic Utilities
- Evaluates exact special values of relative \( L \)-series:
  \[
  L(E/F, 1-j)
  \]
- Uses **PARI/GP via SageMath** for exact arithmetic and special value computations.

---

## File Structure

| File | Description |
| :--- | :--- |
| `global_herm.sage` | Main entry point for global mass calculations and signature handling. |
| `orb_int_herm.sage` | Engine for computing and caching local orbital integrals. |
| `hermi.sage` | Enumeration of Hermitian lattice types and trace matrix precomputation. |
| `jordan.sage` | Jordan decomposition for \( p \)-adic Hermitian and quadratic forms. |
| `quadclass.sage` | Invariants and classification of quadratic forms over \( \mathbb{Q}_p \) and \( \mathbb{Z}/4\mathbb{Z} \). |
| `Lfunc.sage` | Utilities for exact rational evaluation of \( L \)-functions and zeta values. |
| `hermi_density.sage` | Theoretical formulas for local mass terms of Hermitian lattices. |
| `quadlatticeclass.sage` | Classes for handling \( p \)-adic quadratic lattices and their reductions. |
| `deep_extensions_quad_3.sage` | Specialized invariants for quadratic forms over \( \mathbb{F}_3 \). |
| `utils.sage` | Field arithmetic, cyclotomic polynomial factoring, and conjugacy class generation. |

---

## Requirements

- **SageMath**  
  Required for:
  - algebraic number field computations
  - matrix and lattice operations
  - exact arithmetic

- **PARI/GP**  
  Used internally through SageMath for:
  - exact special value calculations
  - \( L \)-function and zeta-related computations

---

## Usage

## 1. Run the Integration Test

The project includes a test script `test_mass_formula.sage` demonstrating the mass computation for \( U(2,2) \) with discriminant \( D_K = -4 \).

```python
load("test_mass_formula.sage")
run_final_integration_test()
```

---

## 2. Manual Mass Computation

You can compute the mass list for a unitary group of dimension \( n \) over a field determined by \( D_K \) as follows.

### Example: \( U(2,2) \) over \( \mathbb{Q}(i) \)

```python
# Example for U(2,2) over Q(i)

n = 4
D_K = -4
primes_list = [2, 3, 5]

# Update the local orbital integral database first
update_orb_int_unitary_db(n, primes_list, D_K)

# Results is a list of [conjugacy_class, mass]
results = mass_list_unitary(n, primes_list, orb_int_unitary_db, D_K)

print(results)
```

---

## Typical Workflow

A standard computation typically proceeds as follows:

1. **Choose the global field / discriminant** \( D_K \)
2. **Specify the rank** \( n \)
3. **Select the relevant local primes**
4. **Precompute local orbital integrals**
5. **Assemble the global mass formula**
6. **Extract the stable conjugacy class masses**

---

## Mathematical Scope

This codebase is intended for computations involving:

- Hermitian lattices over local and global fields
- Stable conjugacy classes of unitary groups
- Local densities and orbital integrals
- Quadratic form classification over \( p \)-adic rings
- Exact arithmetic for special values of \( L \)-functions

---

## References

1. **Olivier Taïbi**,  
   *Dimensions of spaces of level one automorphic forms for split classical groups using the trace formula*,  
   *Ann. Sci. Éc. Norm. Supér.* **50** (2017), no. 2, 269–344.  
   DOI: [10.24033/asens.2323](https://doi.org/10.24033/asens.2323)

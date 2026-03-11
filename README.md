# QUBO Analysis for Molecular Docking

## Table of Contents

- [Overview](#overview)
- [Background: QUBO Formulation](#background-qubo-formulation)
- [Molecular Docking and QUBO Encoding](#molecular-docking-and-qubo-encoding)
- [Solvers Used](#solvers-used)
- [Experimental Results](#experimental-results)
- [Environment Setup](#environment-setup)
- [Code Structure](#code-structure)
- [References](#references)

---

## Overview

This project implements and benchmarks **Quadratic Unconstrained Binary Optimization (QUBO)** models for the purpose of **quantum molecular docking** — a computational drug discovery technique that determines the binding pose and free energy (ΔG_bind) of a small-molecule ligand interacting with a target protein.

Classical docking tools (AutoDock, AutoDock Vina, Glide) rely on heuristic strategies that may sacrifice optimality. This work explores encoding the docking problem as a QUBO instance, enabling the use of quantum annealers and quantum-inspired solvers to search for low-energy binding configurations more systematically.

QUBO is mathematically equivalent to the **Ising model** from statistical physics, positioning it at the intersection of optimization theory and quantum hardware platforms (D-Wave, IBM, Coherent Ising Machines).

---

## Background: QUBO Formulation

A QUBO problem is defined as:

```
minimize  y = xᵀQx
```

where `x` is a binary decision vector and `Q` is a square matrix of constants.

### Key Concepts

- **Diagonal entries** of Q encode linear terms (since xᵢ² = xᵢ for binary variables).
- **Off-diagonal entries** encode quadratic interaction terms between variable pairs.
- **Constraint encoding** — constraints are converted to quadratic penalty terms added to the objective, penalizing infeasible solutions while leaving feasible ones unaffected.
- **Ising equivalence** — via the substitution xᵢ = (1 + sᵢ)/2 with sᵢ ∈ {−1, +1}, any QUBO maps to an Ising Hamiltonian, enabling execution on quantum annealing hardware.

## Molecular Docking and QUBO Encoding

### Problem Setup

Given a ligand with n atoms, molecular docking minimizes the binding free energy ΔG_bind over all possible ligand poses within a docking box D. This is an NP-hard problem.

To encode it as QUBO, the docking box is **discretized into N grid points**, and binary variables `xᵢⱼ ∈ {0, 1}` indicate whether atom aᵢ is matched to grid point gⱼ. Constraints (geometric compatibility, one-atom-per-grid-point) are encoded as Lagrange penalty terms.

Two encoding strategies are implemented:

### Grid Point Matching (GPM)
- Grid points generated with a **2 Å gap** inside the docking box.
- Interaction weights defined using **van der Waals energies** from AutoGrid (AutoDockFR).
- Matches with positive interaction energy are pruned to reduce qubit consumption.

### Feature Atom Matching (FAM)
- Grid points with a **1 Å gap**, coarse-grained into 3 feature atom types: **Hydrophobic (C)**, **H-bond Donor (HD)**, and **H-bond Acceptor (OA)** using the AutoSite algorithm.
- Interaction weights based on **Pauling electronegativity** differences between ligand atoms and feature atoms.
- Docking is split into **pose sampling** (NP-hard, handled by QUBO) and **pose scoring** (classical).

### Encoding Limitations

| Encoding | Max Ligand Atoms | Approx. Peptide Length |
|----------|-----------------|------------------------|
| GPM      | ~156 atoms      | ~15 residues           |
| FAM      | ~305 atoms      | ~30 residues           |

Both methods are limited by **quadratic qubit scaling** with respect to the number of ligand atoms and are not suitable for protein-protein docking.

---

## Solvers Used

| Solver | Type | Backend | Notes |
|--------|------|---------|-------|
| **SCIP** | Exact (MIP/MINLP) | Pyomo | Highly flexible, academic solver with branch-cut-and-price |
| **CPLEX** | Exact (MIP/MIQP) | Pyomo | Commercial solver; branch-and-bound with advanced presolve |
| **Hercules** | Heuristic (QUBO) | Rust library | Lightweight QUBO toolkit; fast but approximate |
| **Qubolite** | Exact (brute-force) | Python/NumPy | Suitable for small/moderate instances; ground-truth reference |

All exact solvers were interfaced through a **Pyomo**-based modelling framework, with QUBO instances passed as binary quadratic programs.

---

## Experimental Results

### Benchmark QUBO Matrices (Q1–Q4)

Tested on matrices of size 4×4, 5×5, and 8×8:

- **SCIP, CPLEX, and Qubolite** consistently converged to identical optimal objective values across all instances.
- **Degeneracy confirmed**: Instances Q2 (5×5) and Q4 (8×8) showed **Jaccard Similarity = 0.00** between SCIP and CPLEX bitstrings despite matching objective values, empirically confirming multiple global optima.
- **Hercules** converged quickly but settled at approximate solutions above the global minimum.

| Matrix | Optimal Value | Notes |
|--------|--------------|-------|
| Q1 (4×4) | -11.0 | Unique global optimum |
| Q2 (5×5) | -5.0  | Degenerate — multiple optima |
| Q3 (8×8) | -12.0 | Unique global optimum |
| Q4 (8×8) | -6889.0 | Highly degenerate- multiple optima |

### Molecular Docking QUBO (Protein: 1Y6R, FAM/GPM, 675 variables)

| Solver | Best Objective | Active Variables | Runtime |
|--------|---------------|-----------------|---------|
| SCIP   | **-8.57**     | 18              | >2 hours |
| CPLEX  | **-8.57**     | 18              | >2 hours |
| Hercules (B&B) | -3.235 | 17         | ~10 seconds |

- SCIP and CPLEX reached the same global minimum (-8.57) with **83.3% structural overlap** (15 common active bits), confirming a stable solution backbone with minor degenerate variations.
- Hercules showed **0% structural overlap** with exact solvers, indicating topological divergence into an unrelated local minimum.
- Key finding: **objective value alone is insufficient** to evaluate QUBO solutions — bitstring-level structural analysis is essential, especially in docking where bit positions represent spatial configurations.

---

## Environment Setup

Two install paths are available depending on which solvers you need.
SCIP and CPLEX are **conda-only** and cannot be installed via pip — Path A is required if you need them.

| Path | Solvers available | Requires |
|------|------------------|----------|
| **A — Conda** (recommended) | SCIP, CPLEX, Qubolite, Pyomo | Conda |
| **B — Pip** | Qubolite, Pyomo only | pip / any venv |

---

### Path A — Conda (full install, all solvers)

**Prerequisites:** [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

```bash
# 1. Create and activate the environment
conda create -n qubo-env python=3.11
conda activate qubo-env

# 2. Install conda-only solvers
conda install -c conda-forge scip
conda install -c ibmdecisionoptimization cplex

# 3. Install this repo — all remaining dependencies are pulled in automatically
pip install git+https://github.com/allu0786ansari/Qubo_Analysis.git

# 4. Register the Jupyter kernel
qubo-setup-kernel

# 5. Launch Jupyter
jupyter notebook
```

Select **Python (qubo-env)** as the kernel, then open any notebook under `1_Code/`.

---

### Path B — Pip only (Qubolite + Pyomo, no SCIP/CPLEX)

```bash
# 1. Create and activate a virtual environment
python -m venv qubo-env

# Windows
qubo-env\Scripts\activate
# Mac / Linux
source qubo-env/bin/activate

# 2. Install this repo — all dependencies are pulled in automatically
pip install git+https://github.com/allu0786ansari/Qubo_Analysis.git

# 3. Register the Jupyter kernel
qubo-setup-kernel

# 4. Launch Jupyter
jupyter notebook
```

> ⚠️ **SCIP and CPLEX are not available on this path.** Notebooks using those solvers will raise an import error. Use Path A if you need them.

---

### Hercules

## Code Structure

```
qubo-molecular-docking/
├── pyproject.toml              # Package metadata and all pip dependencies
├── README.md
├── qubo_docking/               # Installable Python package
│   ├── __init__.py
│   └── cli.py                  # qubo-setup-kernel entry point
└── 1_Code/                     # All implementation notebooks and data
    ├── 1y6r_Qubo/
    │   ├── QUBO_1y6R_Pyomo_SCIP.ipynb
    │   └── Qubo_1y6r_Pyomo-CPLIX.ipynb
    ├── Hercules/
    │   ├── Data/
    │   │   └── QUBO_Matrix.txt
    │   └── Qubo_Matrix.ipynb
    ├── Pyomo/
    │   ├── Data/
    │   │   └── matrix.txt
    │   └── QUBO_Matrix.ipynb
    ├── Qubolite/
    │   ├── Data/
    │   │   └── matrix.txt
    │   └── Qubolite.ipynb
    └── QUBO_Tutorial.ipynb     # 1Y6R molecular docking QUBO (GPM/FAM encoding)
```

---

## References
1. QDock — [JinyinZha/QDock on GitHub](https://github.com/JinyinZha/QDock/tree/main)
2. Glover et al. — [A Tutorial on Formulating and Using QUBO Models (arXiv:1811.11538)](https://arxiv.org/abs/1811.11538)
3. J. Zha et al. — "Encoding Molecular Docking for Quantum Computers," *J. Chem. Theory Comput.*, vol. 19, no. 24, pp. 9018–9024, Dec. 2023. [DOI: 10.1021/acs.jctc.3c00943](https://doi.org/10.1021/acs.jctc.3c00943)
4. IBM — [CPLEX Optimization Studio Documentation](https://www.ibm.com/docs/en/icos/22.1.2)
5. SCIP — [SCIP Doxygen Documentation](https://www.scipopt.org/doc/html/)
6. PySCIPOpt — [PySCIPOpt Documentation](https://pyscipopt.readthedocs.io/en/latest/)
7. Pyomo — [Pyomo Documentation 6.9.5](https://pyomo.readthedocs.io/en/stable/)
8. D-Wave — [Ocean SDK Installation](https://docs.dwavequantum.com/en/latest/ocean/install.html)
9. D-Wave — Problem Formulation Guide, 2022. [www.dwavesys.com](https://www.dwavesys.com)
10. Hercules — [DKenefake/hercules on GitHub](https://github.com/DKenefake/hercules)

# QUBO Analysis for Molecular Docking

## Table of Contents

- [Overview](#overview)
- [Background: QUBO Formulation](#background-qubo-formulation)
- [Molecular Docking and QUBO Encoding](#molecular-docking-and-qubo-encoding)
- [Solvers Used](#solvers-used)
- [Experimental Results](#experimental-results)
- [Environment Setup](#environment-setup)
  - [System Dependencies](#system-dependencies)
  - [AutoDockFR (ADFR Suite)](#autodockfr-adfr-suite)
  - [Python Environment (Pyomo / SCIP / CPLEX / Qubolite)](#python-environment-pyomo--scip--cplex--qubolite)
  - [SCIP and CPLEX — Installation Notes](#scip-and-cplex--installation-notes)
  - [Hercules (Docker)](#hercules-docker)
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

---

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
| Q4 (8×8) | -6889.0 | Highly degenerate — multiple optima |

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

### System Dependencies

Install essential build tools and headers before setting up any Python environment:

```bash
sudo apt update
sudo apt install build-essential python3-pip python3-venv

# Python development headers
sudo apt install python3-dev

# System libraries required by OpenBabel (used by ADFR)
sudo apt install -y libsm6 libxext6 libxrender1
```

---

### AutoDockFR (ADFR Suite)

AutoDockFR is required for grid generation (GPM encoding) and feature atom computation (FAM encoding).

#### 1. Download the ADFR Suite

Download the Linux package `ADFRsuite_x86_64Linux_1.0.tar.gz` from one of the following:

- **Official site:** https://ccsb.scripps.edu/adfr/downloads/
- **Google Drive mirror:** https://drive.google.com/drive/folders/14jfy3ivKrlKVk2FhZGAioCFjMWMiRMJe?usp=sharing

#### 2. Extract and Install

```bash
tar -xvzf ADFRsuite_x86_64Linux_1.0.tar.gz
cd ~/ADFRsuite_x86_64Linux_1.0
bash install.sh
```

> **Alternative if `bash install.sh` fails:**
> ```bash
> chmod -R 755 ADFRsuite_x86_64Linux_1.0
> ./install.sh
> ```

#### 3. Add ADFR Suite to System PATH

Open your shell config file:

```bash
nano ~/.bashrc
```

Scroll to the bottom and add the following line (replace `<your_username>` with your actual username):

```bash
export PATH="$PATH:/home/<your_username>/ADFRsuite_x86_64Linux_1.0/bin"
```

Save and reload:

```bash
# Save: Ctrl+O → Enter → Ctrl+X
source ~/.bashrc
```

---

### IBM CPLEX Installation Guide (Linux)

---

#### Step 1: Download the IBM CPLEX Installer

Download the installer binary from one of the following:

- **Google Drive:** https://drive.google.com/drive/folders/14jfy3ivKrlKVk2FhZGAioCFjMWMiRMJe?usp=sharing


The downloaded file will be named: `cplex_studio2211.linux_x86_64.bin`

---

#### Step 2: Navigate to the Installer Directory

```bash
cd <directory_where_installer_is_saved>
e.g: cd downloads
```

---

#### Step 3: Make the Installer Executable

```bash
chmod +x cplex_studio2211.linux_x86_64.bin
```

---

#### Step 4: Run the Installer

```bash
./cplex_studio2211.linux_x86_64.bin
```

During installation:
- Accept the license agreement.
- Choose the installation directory. The recommended path is:
  ```
  /home/<your_username>/ibm/ILOG/CPLEX_Studio2211

  e.g: /home/allu786ansari/ibm/ILOG/CPLEX_Studio2211

  ```
- Complete the installation process.

---

#### Step 5: Set Up Environment Variables

Open your shell config file:

```bash
nano ~/.bashrc
```

Scroll to the bottom and add the following lines (replace `<your_username>` with your actual username):

```bash
export CPLEX_HOME=$HOME/ibm/ILOG/CPLEX_Studio2211/cplex
export PATH=$PATH:$CPLEX_HOME/bin/x86-64_linux
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CPLEX_HOME/bin/x86-64_linux
```

Save and reload:

```bash
# Save: Ctrl+O → Enter → Ctrl+X
source ~/.bashrc
```

---

#### Step 6: Install the CPLEX Python API

> **Note:** This step requires **Python ≤ 3.10**. Activate a compatible virtual environment before proceeding.

At the end of the CPLEX installation, the installer prints a command to install the Python API. Run it in your terminal — it will look similar to:

```bash
python /home/<your_username>/ibm/ILOG/CPLEX_Studio2211/python/setup.py install

e.g: python /home/allu786ansari/ibm/ILOG/CPLEX_Studio2211/python/setup.py install
```

Replace `<your_username>` with your actual username and adjust the path if you chose a custom installation directory.

---

#### Step 7: Verify the Installation

```bash
python3 -c "import cplex; print(cplex.__version__)"
```

If a version number is printed, the installation was successful.

---

### Python Environment (Pyomo / SCIP / CPLEX / Qubolite)

#### 1. Create and activate the conda environment

```bash
conda create -n dockvenv python=3.10 -y
conda activate dockvenv
```

#### 2. Install the package via pip + git

```bash
git clone https://github.com/allu0786ansari/Qubo_Analysis.git
pip install git+https://github.com/allu0786ansari/Qubo_Analysis.git
```

#### 3. Set up the Jupyter kernel

```bash
python -m ipykernel install --user --name dockvenv --display-name "Python (dockvenv)"
```

#### 4. Install SCIP

```bash
conda install -c conda-forge scip
```

#### 5. Install CPLEX Python API (within conda env)

With `dockvenv` activated, run the setup command printed at the end of the CPLEX installation:

```bash
python /home/<your_username>/ibm/ILOG/CPLEX_Studio2211/python/setup.py install

e.g: python /home/allu786ansari/ibm/ILOG/CPLEX_Studio2211/python/setup.py install
```

#### 6. Launch Jupyter Notebook

```bash
jupyter notebook
```

Click the link printed in the terminal to open in your browser.

#### 7. Select the kernel

In the Jupyter interface, select **Python (dockvenv)** as your kernel.

#### 8. Open notebook files to run

Navigate to the `1_Code/` directory and open the relevant notebooks.

---

### SCIP and CPLEX — Installation Notes

Understanding the limitations of each installation method will help you choose the right approach.

#### SCIP

| Method | Access via Pyomo | Limitations |
|--------|-----------------|-------------|
| `pip install pyscipopt` | ❌ No | Full SCIP functionality, but **cannot be used through Pyomo** |
| `conda install -c conda-forge scip` | ✅ Yes | **Recommended** — full access, no variable limits, Pyomo-compatible |

> **Summary:** SCIP can be installed via pip (`pyscipopt`) without any size limitations, but to use it through Pyomo you must install it via conda.

#### CPLEX

| Method | Access via Pyomo | Limitations |
|--------|-----------------|-------------|
| `pip install cplex` | ✅ Yes | **Community edition only** — limited to 1,000 variables and 1,000 constraints |
| IBM CPLEX Optimization Studio (full) | ✅ Yes | No limits — requires manual installation and PATH configuration |

To use the **full CPLEX** (no variable limit), install IBM CPLEX Optimization Studio manually using the steps in the [IBM CPLEX Installation Guide](#ibm-cplex-installation-guide-linux) above, then install the Python API inside your environment:

```bash
python /home/<your_username>/ibm/ILOG/CPLEX_Studio2211/python/setup.py install
```

> **Note:** Python ≤ 3.10 is required for the CPLEX Python API. The recommended `dockvenv` conda environment uses Python 3.10 and is compatible.

> **Summary:** `pip install cplex` gives only the community edition (1,000-variable limit). For the full solver, manual setup via IBM CPLEX Optimization Studio is required.

---

### Hercules (Docker)

Hercules runs inside a Docker container. Follow these steps to set it up.

#### 1. Install Docker

Follow the official Docker installation guide for your OS: https://docs.docker.com/engine/install/

#### 2. Add your user to the Docker group

```bash
sudo usermod -aG docker $USER
newgrp docker
```

#### 3. Test Docker

```bash
docker run hello-world
```

#### 4. Pull and run the Hercules container

```bash
docker run --platform linux/amd64 -it \
  -p 8888:8888 \
  -v /home/allu786ansari/1_Code/Hercules:/workspace \
  dkenefake/hercules \
  bash
```

> **Note:** Adjust the `-v` volume path to match your local `Hercules/` directory.

#### 5. Inside the container — verify and navigate

```bash
ls /workspace
cd /workspace
```

#### 6. Install Jupyter and dependencies

```bash
pip install notebook ipykernel numpy scipy
```

#### 7. Launch Jupyter Notebook

```bash
jupyter notebook --ip=0.0.0.0 --no-browser --allow-root
```

#### 8. Open in browser

Click the link printed in the terminal, for example:

```
http://127.0.0.1:8888/tree?token=<your_token>
```

---

## Code Structure

All implementation notebooks and data are located in the `1_Code/` directory, organized by solver and problem type:

```
1_Code/                         # All implementation notebooks and data
├── 1y6r_Qubo/                  
│   ├── QUBO_1y6R_Pyomo_SCIP.ipynb
│   └── Qubo_1y6r_Pyomo-CPLIX.ipynb
├── Hercules/                   # Hercules QUBO solver implementations
│   ├── Data/
│   │   └── QUBO_Matrix.txt
│   └── Qubo_Matrix.ipynb
├── qubo_analysis               # Package source
│   └── __init__.py
├── Pyomo/                      # Pyomo + SCIP/CPLEX implementations
│   ├── Data/
│   │   └── matrix.txt
│   └── QUBO_Matrix.ipynb
├── Qubolite/                   # Qubolite solver implementations
│   ├── Data/
│   │   └── matrix.txt
│   └── Qubolite.ipynb
├── qubo_analysis/                   # Qubolite solver implementations
│   ├── __init__.py

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

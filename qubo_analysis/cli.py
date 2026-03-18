"""
CLI entry points for qubo-molecular-docking.
"""

import subprocess
import sys


def setup_kernel():
    """Register the conda environment as a Jupyter kernel."""
    cmd = [
        sys.executable, "-m", "ipykernel", "install",
        "--user",
        "--name", "qubo-env",
        "--display-name", "Python (qubo-env)",
    ]
    print("Registering Jupyter kernel 'Python (qubo-env)' ...")
    result = subprocess.run(cmd, check=True)
    if result.returncode == 0:
        print("Kernel registered successfully.")
        print("Launch Jupyter with:  jupyter notebook")
        print("Then select kernel:   Python (qubo-env)")


if __name__ == "__main__":
    setup_kernel()
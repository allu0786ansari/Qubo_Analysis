.PHONY: install kernel notebook clean

## Full environment setup (recommended)
install:
	conda env create -f environment.yml
	@echo ""
	@echo "========================================"
	@echo "  Environment created: qubo-env"
	@echo "  Activate with:  conda activate qubo-env"
	@echo "  Then run:       make kernel"
	@echo "========================================"

## Update an already-existing environment
update:
	conda env update -f environment.yml --prune

## Register Jupyter kernel (run after activating qubo-env)
kernel:
	qubo-setup-kernel

## Launch Jupyter Notebook
notebook:
	jupyter notebook

## Remove the conda environment entirely
clean:
	conda env remove -n qubo-env
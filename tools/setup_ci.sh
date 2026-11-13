#!/bin/bash -e

# Set up an environment to run CI tests, e.g. with GitHub Actions or Travis

eval "$(conda shell.bash hook)"
conda activate test
conda install --yes pip scipy matplotlib pandas pyrmsd imp-nightly cmake gnuplot

pip install pytest-cov coverage pytest-flake8 'pytest<9'

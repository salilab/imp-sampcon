#!/bin/bash -e

# Set up an environment to run CI tests, e.g. with GitHub Actions or Travis

if [ $# -ne 1 ]; then
  echo "Usage: $0 python_version"
  exit 1
fi

python_version=$1

# Don't install conda's gnuplot on Python 2.7; it conflicts with the IMP package
if [ "${python_version}" != "2.7" ]; then
  gnuplot="gnuplot"
  pip="pip"
else
  pip="pip<=19.3.1"
fi

conda config --remove channels defaults  # get conda-forge, not main, packages
conda create --yes -q -n python${python_version} -c salilab -c conda-forge python=${python_version} ${pip} scipy matplotlib pandas pyrmsd imp-nightly cmake ${gnuplot}
eval "$(conda shell.bash hook)"
conda activate python${python_version}

if [ ${python_version} = "2.7" ]; then
  # pytest-flake8 1.1.0 tries to import contextlib.redirect_stdout, which
  # isn't present in Python 2
  pip install pytest-cov 'coverage<=6.2' 'pytest-flake8<1.1'
else
  pip install pytest-cov coverage pytest-flake8
fi

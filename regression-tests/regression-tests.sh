#!/bin/bash

mkdir -p regression-tests/reference_data
git clone https://github.com/irmb/test_data regression-tests/reference_data

python3 -m venv .venv
source .venv/bin/activate

pip install git+https://github.com/soerenPeters/meshio@update-pyproject-version
pip install git+https://gitlab.com/dglaeser/fieldcompare


./regression-tests/driven_cavity_test.sh
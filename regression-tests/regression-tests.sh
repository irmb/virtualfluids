#!/bin/bash

#################################
# VirtualFludis regression tests
#################################


# 1. Cloning the reference data from github
rm -r reference_data && mkdir -p reference_data
git clone https://github.com/irmb/test_data reference_data

# 2. set up the python environnement
#    by cloning our meshio patch and fieldcompare into a venv
python3 -m venv .venv
source .venv/bin/activate
pip install fieldcompare

# 3. Running the specific tests
./regression-tests/driven_cavity_uniform_test.sh
./regression-tests/driven_cavity_test.sh



# How to add a new regression test?
# 1. setup the specfic simulation and run it to create reference data.
# 2. fork https://github.com/irmb/test_data and create a pull request containing the reference data.
# 3. copy ./regression-tests/driven_cavity_test.sh and adjust the file accordingly to the new test scenario.
# 4. execute this file from here accordingly to #3.
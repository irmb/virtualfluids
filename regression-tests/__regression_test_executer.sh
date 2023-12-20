# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

download_reference_data () {
    rm -rf reference_data && mkdir -p reference_data
    git clone --depth 1 --filter=blob:none --sparse https://github.com/irmb/test_data reference_data
    cd reference_data
    git sparse-checkout add $1
    cd ..
}

# run regression test - arguments:
# 1. REFERENCE_DATA_DIR - to download the reference data and compare against
# 2. CMAKE_FLAGS - cmake flags for the build of VirtualFluids
# 3. APPLICATION - the application to be executed
# 4. RESULT_DATA_DIR - the path to the produced data to be compared
run_regression_test () {
    download_reference_data $1

    rm -rf build && mkdir -p build
    cmake -B build $2
    cmake --build build --parallel 8

    # execute the application
    $3

    # execute fieldcompare (A more comprehensive manual can be found here https://gitlab.com/dglaeser/fieldcompare)
    fieldcompare dir $4 reference_data/$1 --include-files "*.vtu" $5
}

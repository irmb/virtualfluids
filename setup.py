import sys
from pathlib import Path

import setuptools
import skbuild

"""
Install python wrapper of Virtual Fluids
install via python:
    python setup.py install
    set CMAKE Flags via -DBUILD_VF_GPU:BOOL=1
    CMAKE flags have to be separated by -- 
    example: python setup.py install -- VBUILD_VF_CPU:BOOL=ON
or install via pip:
    pip install .
    for pip>21:
        set CMAKE Flags via --config-settings "-DBUILD_VF_GPU=ON"
        example: pip install . --config-settings="-DBUILD_VF_GPU=ON"
        each option has to be passed in individually i.e --config-settings="-DOPT1=ON" --config-settings="-DOPT2=OFF"
    for pip <21:
        set CMAKE Flags via --global-option ="-DBUILD_VF_GPU=ON"
        example: pip install . --global-option="-DBUILD_VF_GPU=ON"
"""

package_name = "pyfluids"
target = "python_bindings"
src_dir = "pythonbindings"

# hack to get config-args for installation with pip>21
cmake_args = []
if("config_args" in locals()):
    cmake_args.extend([f"{k}={v}" for k,v in locals()["config_args"].items()])

cmake_args += [
        f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
        "-DBUILD_VF_PYTHON_BINDINGS=ON",
        "-DBUILD_SHARED_LIBS=OFF",
        "-DBUILD_VF_DOUBLE_ACCURACY=OFF",
        "-DBUILD_VF_UNIT_TESTS:BOOL=OFF",
        "-DBUILD_WARNINGS_AS_ERRORS=OFF",
    ]

skbuild.setup(
    name=package_name,
    packages=[package_name],
    package_dir={"": src_dir},
    cmake_args = cmake_args,
    cmake_install_target=target,
    include_package_data=True,
)

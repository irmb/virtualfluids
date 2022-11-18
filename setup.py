import sys
from pathlib import Path

import skbuild

"""
Install python wrapper of virtual fluids
install via python:
    python setup.py install build_ext
    set CMAKE Flags via -DBUILD_VF_GPU:BOOL=1
    CMAKE flags have to be separated by -- 
    example: python setup.py install build_ext -- VBUILD_VF_CPU:BOOL=ON
or install via pip:
    pip install -e .
    set CMAKE Flags via --config-settings -DBUILD_VF_GPU=1
"""

init_py = "from .bindings import *"
top_dir = Path(__file__).parent

#create __init__.py
pyfluids_dir =  top_dir / "pythonbindings/pyfluids"
pyfluids_dir.mkdir(exist_ok=True)
init_file = pyfluids_dir / "__init__.py"
init_file.write_text(init_py)
(pyfluids_dir / "__init__.py").touch()

# create build_dir
name_of_build_dir = "build"

build_dir = (top_dir/name_of_build_dir).mkdir(exist_ok=True)
extra_args = []
if("cmake_args" in locals()):
    extra_args.extend([f"{k}={v}" for k,v in locals()["cmake_args"].items()])

cmake_args = [
        f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
        "-DBUILD_VF_PYTHON_BINDINGS=ON",
        "-DBUILD_SHARED_LIBS=OFF",
        "-DBUILD_VF_DOUBLE_ACCURACY=OFF",
        "-DBUILD_VF_UNIT_TESTS:BOOL=OFF",
        "-DBUILD_WARNINGS_AS_ERRORS=OFF",
    ] + extra_args

skbuild.setup(
    name="pyfluids",
    packages=["pyfluids"],
    package_dir={"": "pythonbindings"},
    cmake_args = cmake_args,
    cmake_install_dir=str(build_dir),
    cmake_install_target="python_bindings"
)

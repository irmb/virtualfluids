import sys
from pathlib import Path

import setuptools
import skbuild

"""
Install python wrapper of virtual fluids
install via python:
    python setup.py develop
    set CMAKE Flags via -DBUILD_VF_GPU:BOOL=1
    CMAKE flags have to be separated by -- 
    example: python setup.py develop -- VBUILD_VF_CPU:BOOL=ON
or install via pip:
    pip install -e .
    for pip>21:
        set CMAKE Flags via --config-settings "-DBUILD_VF_GPU=ON"
        example: pip install -e . --config-settings="-DBUILD_VF_GPU=ON"
    for pip <21:
        set CMAKE Flags via --global-option ="-DBUILD_VF_GPU=ON"
        example: pip install -e . --global-option="-DBUILD_VF_GPU=ON"
"""

init_py = "from .bindings import *"
top_dir = Path(__file__).parent

#create __init__.py
pyfluids_dir =  top_dir / "pythonbindings/pyfluids"
pyfluids_dir.mkdir(exist_ok=True)
init_file = pyfluids_dir / "__init__.py"
init_file.write_text(init_py)
(pyfluids_dir / "__init__.py").touch()

# # create build_dir
# name_of_build_dir = "build"

# build_dir = top_dir/name_of_build_dir
# build_dir.mkdir(exist_ok=True)


target = "python_bindings"

config_args = []
if("cmake_args" in locals()):
    config_args.extend([f"{k}={v}" for k,v in locals()["cmake_args"].items()])

    
cmake_args = [
        f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
        "-DBUILD_VF_PYTHON_BINDINGS=ON",
        "-DBUILD_SHARED_LIBS=OFF",
        "-DBUILD_VF_DOUBLE_ACCURACY=OFF",
        "-DBUILD_VF_UNIT_TESTS:BOOL=OFF",
        "-DBUILD_WARNINGS_AS_ERRORS=OFF",
    ] + config_args #+ sys_args

# maker = skbuild.cmaker.CMaker()
# maker.configure(clargs=cmake_args, cmake_install_dir=build_dir)
# maker.make(install_target=target)

skbuild.setup(
    name="pyfluids",
    packages=["pyfluids"],
    package_dir={"": "pythonbindings"},
    cmake_args = cmake_args,
    cmake_install_target=target,
    # package_data={"pyfluids": ["bindings.*"]},
    include_package_data=True,
    # cmake_install_dir="pythonbindings"
    # data_files=[("python_bindings/pyfluids", ["bindings.*"])]
)
# setuptools.setup(
#      name="pyfluids",
#     packages=["pyfluids"],
#     package_dir={"": "pythonbindings"},
# )

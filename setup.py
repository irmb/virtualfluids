"""
Install python wrapper of Virtual Fluids
install via python:
    python setup.py install
    set CMAKE Flags via -DBUILD_VF_GPU:BOOL=ON
    CMAKE flags have to be separated by -- 
    example: python setup.py install -- -DBUILD_VF_CPU:BOOL=ON
or install via pip:
    pip install .
    for pip>21:
        set CMAKE Flags via --config-settings "-DBUILD_VF_GPU=ON"
        example: pip install . --config-settings="-DBUILD_VF_GPU=ON"
        each option has to be passed in individually i.e
        --config-settings="-DOPT1=ON" --config-settings="-DOPT2=OFF"
    for pip <21:
        set CMAKE Flags via --global-option ="-DBUILD_VF_GPU=ON"
        example: pip install . --global-option="-DBUILD_VF_GPU=ON"
"""
import sys
from pathlib import Path

from setuptools import find_packages
import skbuild

package_name = "pyfluids"
target = "python_bindings"
src_dir = "pythonbindings"
stubs_package = package_name+"-stubs"
stub_dir = Path(src_dir)/stubs_package


def find_stub_subpackages(stub_dir: Path):
    return [str(d.parent.relative_to(stub_dir.parent)) for d in stub_dir.rglob("__init__.pyi")]


def find_stub_files(dir: Path):
    return [str(f.relative_to(dir)) for f in dir.rglob("*.pyi")]


# hack to get config-args for installation with pip>21
cmake_args = []
if "config_args" in locals():
    cmake_args.extend([f"{k}={v}" for k, v in locals()["config_args"].items()])

cmake_args += [
        f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
        "-DBUILD_VF_PYTHON_BINDINGS=ON",
        "-DBUILD_SHARED_LIBS=OFF",
        "-DBUILD_VF_DOUBLE_ACCURACY=OFF",
        "-DBUILD_VF_UNIT_TESTS:BOOL=OFF",
        "-DBUILD_WARNINGS_AS_ERRORS=OFF",
    ]
print(find_stub_subpackages(stub_dir))
print(find_packages(where=src_dir))
skbuild.setup(
    name=package_name,
    packages=find_packages()+find_stub_subpackages(stub_dir),
    package_dir={"": src_dir},
    cmake_args=cmake_args,
    cmake_install_target=target,
    package_data={package_name: ["py.typed"], stubs_package: find_stub_files(stub_dir)},
    include_package_data=True,
)

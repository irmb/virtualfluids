import sys
from pathlib import Path
from typing import List

import skbuild

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
        each option has to be passed in individually i.e --config-settings="-DOPT1=ON" --config-settings="-DOPT2=OFF"
    for pip <21:
        set CMAKE Flags via --global-option ="-DBUILD_VF_GPU=ON"
        example: pip install . --global-option="-DBUILD_VF_GPU=ON"
"""

package_name = "pyfluids"
target = "python_bindings"
src_dir = "pythonbindings"
stub_package = package_name+"-stubs"

stub_dir = Path(src_dir)/stub_package


def add_subfiles(dir_path: Path, suffix: str, root_dir: Path) -> List[str]:
    files = []
    for f in dir_path.iterdir():
        if f.is_dir():
            files.extend(add_subfiles(f, suffix, root_dir))
        if f.is_file():
            if f.suffix != suffix:
                continue
            files.append(str(f.relative_to(root_dir)))
    return files

def add_directory(dir_path: Path, suffix: str):
    return add_subfiles(dir_path, suffix, dir_path)

stub_files = add_directory(stub_dir, ".pyi")

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

skbuild.setup(
    name=package_name,
    packages=[package_name, "pymuparser", "pyfluids-stubs"],
    package_dir={"": src_dir},
    cmake_args=cmake_args,
    cmake_install_target=target,
    package_data={  "pyfluids": ["py.typed"],
                    "pyfluids-stubs": stub_files},
    include_package_data=True,
)

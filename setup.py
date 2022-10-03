import inspect
import sys
from pathlib import Path

import cmake_build_extension
import setuptools

"""
Install python wrapper of virtual fluids
install via python setup.py install build_ext
set CMAKE Flags via -DBUILD_VF_GPU:BOOL=1
"""

init_py = inspect.cleandoc(
    """
    import cmake_build_extension
    with cmake_build_extension.build_extension_env():
        from .bindings import *
    """
)

extra_args = []
if("cmake_args" in locals()):
    extra_args.extend([f"{k}={v}" for k,v in locals()["cmake_args"].items()])

setuptools.setup(
    ext_modules=[
        cmake_build_extension.CMakeExtension(
            name="pyfluids",
            install_prefix="pyfluids",
            write_top_level_init=init_py,
            source_dir=str(Path(__file__).parent.absolute()),
            cmake_configure_options = [
                f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
                "-DCALL_FROM_SETUP_PY:BOOL=ON",
                "-DBUILD_VF_PYTHON_BINDINGS=ON",
                "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache",
                "-DCMAKE_CUDA_COMPILER_LAUNCHER=ccache",
                "-DCMAKE_C_COMPILER_LAUNCHER=ccache",
                "-DBUILD_SHARED_LIBS=OFF",
                "-DBUILD_VF_DOUBLE_ACCURACY=OFF",
                "-DBUILD_VF_UNIT_TESTS:BOOL=OFF",
                "-DBUILD_WARNINGS_AS_ERRORS=OFF",
            ] + extra_args,
        )
    ],
    cmdclass=dict(
        build_ext=cmake_build_extension.BuildExtension,
    ),
)

import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools.command.develop import develop
from distutils.version import LooseVersion

# Installation via pip: pip install -e .
# If GPU backend: pip install -e . --install-option="--GPU"

vf_cmake_args = [
    "-DBUILD_VF_PYTHON_BINDINGS=ON",
    "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache",
    "-DCMAKE_CUDA_COMPILER_LAUNCHER=ccache",
    "-DCMAKE_C_COMPILER_LAUNCHER=ccache",
    "-DUSE_MPI=ON",
    "-DBUILD_SHARED_LIBS=OFF",
    "-DBUILD_WARNINGS_AS_ERRORS=OFF"
]

vf_cpu_cmake_args = [
    "-DBUILD_VF_DOUBLE_ACCURACY=ON",
    "-DBUILD_VF_CPU:BOOL=ON",
    "-DBUILD_VF_UNIT_TESTS:BOOL=ON",
    "-DUSE_METIS=ON",
]

vf_gpu_cmake_args = [
    "-DBUILD_VF_DOUBLE_ACCURACY=OFF",
    "-DBUILD_VF_GPU:BOOL=ON",
    "-DBUILD_VF_UNIT_TESTS:BOOL=OFF",
    "-DUSE_METIS=OFF",
]


class CommandMixin:
    user_options = [
        ('GPU', None, 'compile pyfluids with GPU backend'),
    ]

    def initialize_options(self):
        super().initialize_options()
        self.GPU = False

    def finalize_options(self):
        super().finalize_options()

    def run(self):
        global GPU
        GPU = self.GPU
        super().run()

class InstallCommand(CommandMixin, install):
    user_options = getattr(install, 'user_options', []) + CommandMixin.user_options

class DevelopCommand(CommandMixin, develop):
    user_options = getattr(develop, 'user_options', []) + CommandMixin.user_options

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(CommandMixin, build_ext):
    user_options = getattr(build_ext, 'user_options', []) + CommandMixin.user_options

    def run(self):
        super().run()
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        cmake_args.extend(vf_cmake_args)
        cmake_args.extend(vf_gpu_cmake_args if GPU else vf_cpu_cmake_args)

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


setup(
    name='pyfluids',
    version='0.0.1',
    ext_modules=[CMakeExtension('pyfluids')],
    cmdclass={"install": InstallCommand, "develop": DevelopCommand, "build_ext":CMakeBuild},
    zip_safe=False,
)

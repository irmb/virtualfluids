from skbuild import setup

cmake_args = ["-DBUILD_VF_CPU:BOOL=ON", "-DUSE_METIS=ON", "-DUSE_MPI=ON", "-DBUILD_SHARED_LIBS=OFF",
              "-DBUILD_VF_UNIT_TESTS:BOOL=ON"]

setup(
    name="virtualfluids",
    version="0.0.1",
    author="Sven Marcus",
    author_email="sven.marcus@tu-braunschweig.de",
    cmake_args=cmake_args
)

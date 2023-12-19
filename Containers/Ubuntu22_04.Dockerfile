# VirtualFluids Development Image:
# Ubuntu 22.04

FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

# timezone
ARG TZ
ENV TZ="$TZ"

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update &&   \
    apt-get install -y  \
    git wget unzip software-properties-common \
    build-essential g++-12 gcc-12 gdb \
    ccache              \
    ninja-build         \
    openmpi-bin         \
    libopenmpi-dev      \
    libomp-15-dev       \
    clang-15            \
    clang-format-15     \
    clang-tidy-15       \
    clang-tools-15      \
    python3.11          \
    python3-pip         \
    python3.11-dev      \
    cppcheck            \
    # needed for doxygen
    flex bison          \ 
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 \
    && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-15 100 \
    && update-alternatives --install /usr/bin/clang clang /usr/bin/clang-15 100 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 100 \
    && ln -s clang-tidy-15 /usr/bin/clang-tidy \
    && ln -s clang-format-15 /usr/bin/clang-format \
    && wget https://github.com/clangd/clangd/releases/download/16.0.2/clangd-linux-16.0.2.zip && unzip clangd-linux-16.0.2.zip \
    && mv clangd_16.0.2/bin/clangd /usr/bin/clangd-16 && mv clangd_16.0.2/lib/clang/16 /usr/lib/clang/ \
    && update-alternatives --install /usr/bin/clangd clangd /usr/bin/clangd-16 100 \
    && pip3 install      \
        cmake==3.26.3    \
        setuptools       \
        wheel            \
        scikit-build     \
        pyvista          \
        numpy            \
        ansible          \
        'jinja2<3.1'     \
        gcovr==6.0       \
        lizard==1.17.10  \
    && wget https://www.doxygen.nl/files/doxygen-1.9.8.src.tar.gz && tar -xvf doxygen-1.9.8.src.tar.gz \
    && cd doxygen-1.9.8 && mkdir build && cd build && cmake -G "Unix Makefiles" .. && make -j8 && make install
# flex and bison apt packages are needed for doxygen

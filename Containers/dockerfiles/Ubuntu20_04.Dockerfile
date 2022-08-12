# VirtualFluids BuildDependencies:
# Ubuntu 20.04
# general tools: wget, unzip, git
# CMake 3.24.0
# ccache
# gcc 9.3 (default)
# gdb
# openmpi 4.0.3
# openmp
# cuda 11.3.1 as base image
# freeGLUT
# cppcheck
# clang 10.0 (default)
# clangd language server https://clangd.llvm.org/
# python pip3 with modules: setuptools, wheels, scikit-build, pyvista, numpy, ansible, gcovr

FROM nvidia/cuda:11.3.1-devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update &&   \
    apt-get install -y  \
    wget unzip git      \
    build-essential gdb \
    ccache              \
    ninja-build         \
    openmpi-bin=4.0.3-0ubuntu1 \
    libomp-dev          \
    clang               \
    clang-format        \
    clang-tidy          \
    clang-tools         \
    llvm-dev            \
    libclang-dev        \
    python3-pip         \
    freeglut3-dev       \
    cppcheck            \
    clangd-12           \
    && update-alternatives --install /usr/bin/clangd clangd /usr/bin/clangd-12 100 \
    && mkdir -p /usr/local/cmake/ && cd /usr/local/cmake/ \
    && version=3.24 && build=0 \
    && wget https://cmake.org/files/v$version/cmake-$version.$build-linux-x86_64.tar.gz \
    && tar -xzvf cmake-$version.$build-linux-x86_64.tar.gz                              \
    && ln -s /usr/local/cmake/cmake-$version.$build-linux-x86_64/bin/* /usr/local/bin/  \
    && pip3 install      \
        setuptools       \
        wheel            \
        scikit-build     \
        pyvista          \
        numpy            \
        ansible          \
        'jinja2<3.1' gcovr==5.0

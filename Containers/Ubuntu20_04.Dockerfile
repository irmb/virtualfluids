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

# software-properties-common for add-apt-repository

FROM nvidia/cuda:11.3.1-devel-ubuntu20.04

# timezone
ARG TZ
ENV TZ="$TZ"

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update &&   \
    apt-get install -y  \
    wget unzip software-properties-common \
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
    && pip3 install      \
        cmake==3.26.3    \
        setuptools       \
        wheel            \
        scikit-build     \
        pyvista          \
        numpy            \
        ansible          \
        'jinja2<3.1' gcovr==5.0 \
    && apt update && add-apt-repository -y ppa:git-core/ppa && apt update && apt install git -y

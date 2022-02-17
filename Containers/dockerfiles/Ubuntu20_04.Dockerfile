# VirtualFluids BuildDependencies:
# Ubuntu 20.04
# general tools: wget, unzip, git
# CMake 3.20.5
# ccache
# gcc 9.3 (default)
# openmpi 4.0.3
# openmp
# cuda 11.3.1 as base image
# clang 10.0 (default)

FROM nvidia/cuda:11.3.1-devel-ubuntu20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update &&\
    apt-get install -y\
    wget unzip git \
    build-essential gdb\
    ccache \
    ninja-build \
    openmpi-bin=4.0.3-0ubuntu1 \
    libomp-dev \
    clang clang-format clang-tidy clang-tools llvm-dev libclang-dev &&\
    mkdir -p /usr/local/cmake/ && cd /usr/local/cmake/ &&\
    version=3.20 && build=5 &&\
    wget https://cmake.org/files/v$version/cmake-$version.$build-linux-x86_64.tar.gz &&\
    tar -xzvf cmake-$version.$build-linux-x86_64.tar.gz &&\
    ln -s /usr/local/cmake/cmake-$version.$build-linux-x86_64/bin/* /usr/local/bin/
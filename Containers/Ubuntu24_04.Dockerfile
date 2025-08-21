# #######################################################################################
# ____          ____    __    ______     __________   __      __       __        __
# \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
#  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
#   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
#    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
#     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
#      \    \  |    |   ________________________________________________________________
#       \    \ |    |  |  ______________________________________________________________|
#        \    \|    |  |  |         __          __     __     __     ______      _______
#         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
#          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
#           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
#            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
#
#  This file is part of VirtualFluids. VirtualFluids is free software: you can
#  redistribute it and/or modify it under the terms of the GNU General Public
#  License as published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  SPDX-License-Identifier: GPL-3.0-or-later
#  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
#
# #################################################################################
# VirtualFluids Development Image:
# Ubuntu 24.04

FROM nvidia/cuda:12.9.1-devel-ubuntu24.04 AS base

# timezone
ARG TZ
ENV TZ="$TZ"
ENV USERNAME=ubuntu
ENV GID=1000

ENV VIRTUAL_ENV=/home/$USERNAME/venv

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update &&   \
    apt-get install -y  \
    git git-lfs wget unzip software-properties-common \
    build-essential g++-14 gcc-14 gdb \
    ccache              \
    ninja-build         \
    openmpi-bin         \
    libopenmpi-dev      \
    libomp-18-dev       \
    clang-18            \
    clang-format-18     \
    clang-tidy-18       \
    clang-tools-18      \
    clangd-18      \
    python3.12          \
    python3-pip         \
    python3.12-venv     \
    python3.12-dev      \
    cppcheck            \
    # needed for doxygen
    flex bison 

RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 \
    && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-18 100 \
    && update-alternatives --install /usr/bin/clang clang /usr/bin/clang-18 100 \
    && ln -s clang-tidy-18 /usr/bin/clang-tidy \
    && ln -s clang-format-18 /usr/bin/clang-format \
    && ln -s clangd-18 /usr/bin/clangd

RUN pip3 install cmake==3.26.3 --break-system-packages

RUN wget https://www.doxygen.nl/files/doxygen-1.9.8.src.tar.gz && tar -xvf doxygen-1.9.8.src.tar.gz \
    && cd doxygen-1.9.8 && mkdir build && cd build && cmake -G "Unix Makefiles" .. && make -j8 && make install
    
USER $USERNAME:$GID

RUN python3 -m venv $VIRTUAL_ENV \
    && $VIRTUAL_ENV/bin/python3 -m pip install --upgrade pip \
    && $VIRTUAL_ENV/bin/python3 -m pip install      \
        setuptools       \
        wheel            \
        scikit-build     \
        pyvista          \
        numpy            \
        ansible          \
        'jinja2<3.1'     \
        gcovr==6.0       \
        lizard==1.17.10  \
        matplotlib       \
        mpi4py

FROM base AS wifi
RUN $VIRTUAL_ENV/bin/python3 -m pip install      \
    cupy-cuda12x      \
    git+https://source.coderefinery.org/Hkorb/wifi.git
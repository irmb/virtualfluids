# VirtualFluids Development Image:
# Ubuntu 22.04

FROM git.rz.tu-bs.de:4567/irmb/virtualfluids/ubuntu22_04:1.1 as build

# For wiFI https://source.coderefinery.org/Hkorb/wifi.git
RUN pip3 install      \
    pandas            \
    cupy-cuda12x      \
    mpi4py            \
    && git clone --recurse-submodules https://source.coderefinery.org/Hkorb/wifi.git \
    && cd wifi  \
    && git checkout develop \
    && pip3 install -e . >> log_wifi_install.txt

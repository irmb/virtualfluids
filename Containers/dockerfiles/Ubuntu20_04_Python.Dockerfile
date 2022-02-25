FROM irmb/virtualfluids-deps-ubuntu20.04:latest

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update &&    \
    apt-get install -y   \
        python3          \
        python3-venv     \
        python3-pip      \
        uuid-dev         \
        libgpgme-dev     \
        squashfs-tools   \
        libseccomp-dev   \
        pkg-config       \
        cryptsetup-bin   \
        golang           \
        libgl1           \
    && pip3 install      \
        setuptools       \
        wheel            \
        scikit-build     \
        pyvista          \
        numpy            \
        ansible          \
    && export VERSION=3.8.0 \
    && wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz \
    && tar -xzf singularity-ce-${VERSION}.tar.gz

WORKDIR /singularity-ce-3.8.0
RUN ./mconfig && \
    make -C ./builddir && \
    make -C ./builddir install

WORKDIR /
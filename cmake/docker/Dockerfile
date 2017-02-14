FROM debian:jessie

RUN apt-get update && apt-get install -y build-essential checkinstall wget
RUN set -xe; \
    wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz; \
    tar xf cmake-3.2.2.tar.gz; \
    cd cmake-3.2.2; \
    ./configure; \
    make -j $(nproc); \
    checkinstall

ENV SP2_DEPS \
    cmake \
    make \
    automake \
    build-essential \
    libboost-all-dev \
    libopenmpi-dev

ENV PHONOPY_DEPS \
    python-dev \
    python-numpy \
    python-matplotlib \
    python-yaml \
    python-h5py \
    python-pip

RUN set -xe; \
    \
    apt-get update; \
    \
    LAMMPS_DEPS=$(apt-cache depends lammps | \
        grep Depends | sed "s/.*ends:\ //" | tr '\n' ' '); \
    \
    apt-get install -y    \
            $SP2_DEPS     \
            $PHONOPY_DEPS \
            $LAMMPS_DEPS  \
            apt-file      \
            clang         \
        --no-install-recommends; \
    \
    apt-file update; \
    \
    pip install phonopy

COPY ./build_package.sh /usr/local/bin/

CMD ["build_package.sh"]

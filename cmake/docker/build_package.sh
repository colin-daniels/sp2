#!/bin/bash
set -xe
cd /tmp/sp2/build

package_name=sp2
package_ver=1.0-1
package_loc=$(pwd)/dpkg/${package_name}_${package_ver}
mkdir -p ${package_loc}

# make the executable
CC=clang CXX=clang++ cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DSP2_ENABLE_TESTS=OFF  \
    -DSP2_ENABLE_MPI=ON     \
    -DSP2_ENABLE_PHONOPY=ON \
    -DSP2_BUILD_LIB=OFF     \
    -DBUILD_SHARED_LIBS=OFF

make -j $(nproc)

# install to package dir
make DESTDIR=${package_loc} install

# determine runtime dependencies
sp2_libs=$(ldd sp2 | grep '=>' | sed 's@.\+\w => \([\w/].\+\) (0x\w\+)@\1@')
sp2_deps=$(dpkg -S ${sp2_libs} | sed 's/^\([^:]\+\).\+$/\1/' | sort -u \
    | tr '\n' ' ' | sed 's/ \b/, /g')

mkdir -p ${package_loc}/DEBIAN
cat << EOF > ${package_loc}/DEBIAN/control
Package: ${package_name}
Version: ${package_ver}
Section: base
Priority: optional
Architecture: amd64
Depends: ${sp2_deps}
Maintainer: Colin Daniels <colin.r.daniels@gmail.com>
Description: SP2
EOF

dpkg-deb --build ${package_loc}

#!/bin/sh
set -xe

export MSANDERHOME=`pwd`
./configure --conda --openmp

cd src
make -f Makefile.ap install
cd ..

RSYNC=rsync
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" == "1" ]]; then
  RSYNC=${BUILD_PREFIX}/bin/rsync
fi
${RSYNC} -av bin dat lib $PREFIX

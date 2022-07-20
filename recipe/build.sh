#!/bin/sh
set -xe

# cross compiling options for conda
RSYNC=rsync
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" == "1" ]]; then
  RSYNC=${BUILD_PREFIX}/bin/rsync
  CONDA_SUBDIR=osx-64 mamba install -y -p ${BUILD_PREFIX} -c phenix-project amber_phenix
fi

export MSANDERHOME=`pwd`
./configure --conda --openmp

cd src
make -f Makefile.ap install
cd ..

${RSYNC} -av bin dat lib $PREFIX

#!/bin/sh

export MSANDERHOME=`pwd`
./configure --conda --openmp --verbose

cd src
make -f Makefile.ap install
cd ..

rsync -av bin dat lib $PREFIX

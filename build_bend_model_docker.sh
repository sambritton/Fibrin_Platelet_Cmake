#!/bin/bash

set -ie

if [[ -z $1 ]]; then echo -e "USAGE\n\t$0 BUILD_NAME" && exit 256; fi
BUILD_NAME=$1

rm -rf $BUILD_NAME

mkdir -p $BUILD_NAME && cd $BUILD_NAME

cmake ../src/Linux-Debug/ \
-DCMAKE_C_COMPILER=$(which gcc)  \
-DCMAKE_CXX_COMPILER=$(which g++)  \
-DCUDA_INCLUDE_DIRS=/usr/include/cuda/include \
-DCMAKE_CUDA_COMPILER=$(which nvcc) \
-DBUILD_TESTING=OFF \


make -j 12

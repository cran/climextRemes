#!/bin/bash

SOURCE_DIR=$PWD

climextremes="$SOURCE_DIR/noarch"

PKG_DIR=/tmp/cebuild

if [ ! -e "${PKG_DIR}" ]; then
    mkdir ${PKG_DIR}
fi

# build in a Conda environment
conda create --name climextremes_build python=3.7

source activate climextremes_build
conda install -y conda-build anaconda-client

conda update -y conda
conda update -y conda-build

echo conda build $climextremes
# `conda` should be the conda in the environment but it may be a function
${CONDA_PREFIX}/bin/conda build --output-folder=${PKG_DIR} ${climextremes}

## test with:
# conda install /tmp/cebuild/noarch/climextremes-0.2.1rc7-pyr351_0.tar.bz2

#!/bin/bash

SOURCE_DIR=$PWD

climextremes="$SOURCE_DIR/noarch"

PKG_DIR=/tmp/cebuild

if [ ! -e "${PKG_DIR}" ]; then
    mkdir ${PKG_DIR}
fi

# build in a Conda environment
conda create -y --name climextremes_build python=3.8

source activate climextremes_build
conda install -y conda-build anaconda-client conda-verify

conda update -y conda
conda update -y conda-build

echo conda build $climextremes
# `conda` should be the conda in the environment but it may be a function
# as of 2022-07-13, rpy2 v. 2.9.4 still seems to be latest in anaconda channel
# this is limiting us to pandas <=0.25.3, which is becoming ridiculous.
# So use conda-forge to get newer rpy2 at build and install stages
${CONDA_PREFIX}/bin/conda build -c conda-forge --output-folder=${PKG_DIR} ${climextremes}

## test with:
# conda install /tmp/cebuild/noarch/climextremes-0.2.3rc2-pyr41h7b7c402_0.tar.bz2

#!/bin/bash

SOURCE_DIR=$PWD

climextremes="$SOURCE_DIR/noarch"

PKG_DIR=/tmp/cebuild

if [ ! -e "${PKG_DIR}" ]; then
    mkdir ${PKG_DIR}
fi

# build in a Conda environment
# with 3.11 got issues with numpy/python/pandas versions
mamba create -y --name climextremes_build python=3.10

source activate climextremes_build
mamba install -y conda-build anaconda-client conda-verify

mamba update -y conda
mamba update -y conda-build

echo conda build $climextremes
# `conda` should be the conda in the environment but it may be a function
# as of 2022-07-13, rpy2 v. 2.9.4 still seems to be latest in anaconda channel
# this is limiting us to pandas <=0.25.3, which is becoming ridiculous.
# So use conda-forge to get newer rpy2 at build and install stages
${CONDA_PREFIX}/bin/conda build -c conda-forge --output-folder=${PKG_DIR} ${climextremes}

## test with:
# mamba install /tmp/cebuild/noarch/climextremes-0.3.1rc2-pyr43h4f56d60_0.tar.bz2

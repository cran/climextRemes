#!/bin/bash

SOURCE_DIR=$PWD
CAM_DIR=$1

if [[ -z "${CAM_DIR// }"  ]]; then
  CAM_DIR=/tmp
fi

cd $CAM_DIR

unamestr=`uname`
download_prefix='https://repo.continuum.io/miniconda'
download_file=''

climextremes="$SOURCE_DIR/osx"

if [[ "$unamestr" == 'Linux' ]]; then
   download_file='Miniconda3-latest-Linux-x86_64.sh'
   climextremes="$SOURCE_DIR/linux"
elif [[ "$unamestr" == 'Darwin' ]]; then
   download_file='Miniconda3-latest-MacOSX-x86_64.sh'
   climextremes="$SOURCE_DIR/osx"
else
   printf "Unsupported format $unamestr. Aborting...\n"
   exit 0
fi


printf "Setting up environment on $PWD\n"
curl -s $download_prefix/$download_file > $download_file
bash $download_file -f -b -p $CAM_DIR/conda
source $CAM_DIR/conda/bin/activate
source $CAM_DIR/conda/bin/activate

printf "Finished setting up environment\n"
# output to do sanity check..


conda install -y conda-build anaconda-client

conda update -y conda
conda update -y conda-build

echo conda build $climextremes
conda build $climextremes



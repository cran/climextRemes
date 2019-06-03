PYVER=3.7
CEVER=0.2.1rc1
CHANNEL=cascade # paciorek

Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"

## [DONE] Start with Conda + Python 3.7 + fresh install

conda create -y --name ce_test_${PYVER} python=${PYVER}
source activate ce_test_${PYVER}

conda install -y -c ${CHANNEL} climextremes=${CEVER}
# this pulls in r 3.5.1-mro, r-rcpp, r-rcpparmadillo, numpy 1.16.4, pandas 0.24.2, rpy2 2.9.4, tzlocal 1.5.1

# note that R packages are installed in ~/R/...  and NOT in isolated conda env (ugggh)
Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"
python -c "import climextremes"
python ~/Desktop/test.py

## Now try without R climextRemes and R package dependencies

Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"

python ~/Desktop/test.py

## Now try without R climextRemes 

Rscript -e "remove.packages(c('extRemes','distillery','Lmoments'))"

python ~/Desktop/test.py

## using climextRemes 0.2.1

# Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"
# don't remove dependencies as can get into issues with versions in different libraries (I think)
which R
R CMD INSTALL ~/Desktop/climextRemes_0.2.1.tar.gz

python ~/Desktop/test.py

## Now with old R - can't do on Linux

conda create -y --name ce_test_oldR_${PYVER} python=${PYVER}
source activate ce_test_oldR_${PYVER}

conda install -y r==3.4.3  # can't get 3.4.0
conda install -y -c ${CHANNEL} climextremes=${CEVER}

# UnsatisfiableError: The following specifications were found to be in conflict:
#  - climextremes=0.2.1rc18


Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"
python ~/Desktop/test.py

## Now with new rpy2 - this doesn't work as the climextremes install downgrades rpy2

R CMD INSTALL ~/Desktop/climextRemes_0.2.1.tar.gz

conda uninstall -y rpy2
pip install rpy2==3.0.0
conda install -y -c ${CHANNEL} climextremes=${CEVER}

conda uninstall -y rpy2 climextremes
pip install rpy2==3.0.0
pip install --user --index-url https://test.pypi.org/simple/ climextremes==${CEVER}

python ~/Desktop/test.py
# runs but R pkg installation seems to not find cran or cnr repos, so need R pkgs already installed

## [DONE] Python 3.5

PYVER=3.5

conda create -y --name ce_test_${PYVER} python=${PYVER}
source activate ce_test_${PYVER}

conda install -y -c ${CHANNEL} climextremes=${CEVER}
# this pulls in r 3.5.1-mro, r-rcpp, r-rcpparmadillo, numpy 1.15.2, pandas 0.23.4, rpy2 2.9.4, tzlocal 1.5.1

Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"
python -c "import climextremes"
# got one failure here with Lmoments not being available to extRemes - I think because it was using my own ~/R library by default
python ~/Desktop/test.py


## Python 2.7

PYVER=2.7
CEVER=0.2.1rc17

conda create -y --name ce_test_${PYVER} python=${PYVER}
source activate ce_test_${PYVER}

conda install -y -c ${CHANNEL} climextremes=${CEVER}
# pulls in r 3.5.1mro, numpy 1.16.4, pandas 0.24.2 r-rcpp, r-rcpparmadillo, rpy2 2.8.6, tzlocal 1.5.1 
Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"

python -c "import climextremes"

python ~/Desktop/test.py

## pip install climextremes outside conda
## this won't pull in R; assume user has R installed

pip install --user tzlocal
pip install --user --index-url https://test.pypi.org/simple/ climextremes==${CEVER}

python ~/Desktop/test.py

## now use pip + rpy2 >= 3.0.0

Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"
pip install --user rpy2==3.0.0
python ~/Desktop/test.py  # this installed R  packages (with some extra output) but of course needs 0.2.1 so or fails with NULLType error
R CMD INSTALL ~/Desktop/climextRemes_0.2.1.tar.gz
python ~/Desktop/test.py

## Mac with conda

conda create -y --name ce_test_mac_${PYVER} python=${PYVER}
source activate ce_test_mac_${PYVER}

conda install -y libpng=1.6.36  # unsafe path in 1.6.37
conda install -y -c ${CHANNEL} climextremes=${CEVER}
# this pulls in r 3.5.1, r-rcpp, r-rcpparmadillo, numpy 1.16.4, pandas 0.24.2, rpy2 2.9.4, tzlocal 1.5.1
# downgrades python 3.7.3 to python 3.7.0 for some reason
# SafetyError: The package for r-base located at /accounts/gen/vis/paciorek/.conda/pkgs/r-base-3.5.1-h539fb6c_1 (but seems to install anyway)

# note that R packages are installed in ~/R/...  and NOT in isolated conda env (ugggh)
Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"

# pkg_resources.DistributionNotFound: The 'cffi>=1.0.0' distribution was not found and is required by rpy2 
conda install -y cffi
# that might be because it was actually trying to use rpy2 >=3.0.0

python ~/Desktop/test.py

## Mac with conda (laptop)

conda create -y --name ce_test_mac_${PYVER} python=${PYVER}
source activate ce_test_mac_${PYVER}

conda install -y -c ${CHANNEL} climextremes=${CEVER}
# this pulls in r 3.5.1, r-rcpp, r-rcpparmadillo, numpy 1.16.4, pandas 0.24.2, rpy2 2.9.4, tzlocal 1.5.1
# same safetyError as above

python ~/Desktop/test.py

## Mac with pip and R 3.6 (laptop)

pip install rpy2  # gets >= 3.0.0
# install 0.2.1 climextremes therefore from tarball
pip install --user --index-url https://test.pypi.org/simple/ climextremes==${CEVER}
# pip install --user climextremes==${CEVER}

python ~/Desktop/test.py

## Mac with pip and old R (laptop)

# see above

## Windows on GCP

# follow GCP windows instructions: https://cloud.google.com/compute/docs/quickstart-windows
# turn off enhanced security in server manager
# install miniconda
# start powershell Anaconda prompt

# then create,activate env and install climextremes


                                                                          

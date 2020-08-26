PYVER=3.8
CEVER=0.2.2rc1

# variations:
# Python 3.7, 3.8
# R 3.6, 4.0
# rpy2 2.9.4, 3.2.5, 3.3.5
# conda, pip, conda+conda-R
# w/ R climextRemes missing, 0.2.1, 0.2.2
# w/ and w/o extRemes, distillery, Lmoments
# w/ and w/o rpy2, numpy, pandas


### test PyPI: do on Linux and Linux under 3.8 (don't have virtualenv readily accessible on Mac)

Rscript -e "remove.packages(c('extRemes','climextRemes','distillery','Lmoments'))"
R CMD INSTALL climextRemes_${CEVER}.tar.gz  # otherwise 'import climextremes' will pull in old CRAN version

VE=~/Desktop/ce_${CEVER}_${PYVER}
virtualenv ${VE}
source ${VE}/bin/activate

## this gets pandas >= 1.0.0 and rpy2 >= 3.3.0
pip install rpy2 numpy pandas tzlocal # needed because not in the index-url, presumably
pip install --index-url https://test.pypi.org/simple/ climextremes==${CEVER}

## now test:

python ~/Desktop/test.py

## try alternative rpy2 versions and run the fit_gev_examples.py
pip install rpy2==3.2.5 pandas==1.1.0
pip install rpy2==2.9.4 pandas==0.25.3

### Conda from paciorek channel: do on Linux and Mac

CHANNEL=paciorek

conda create -y --name ce_${CEVER}_${PYVER} python=${PYVER}
source activate ce_${CEVER}_${PYVER}

conda install -y -c ${CHANNEL} climextremes=${CEVER}

## this will use old R climextRemes
python ~/Desktop/test.py

## now try with updated climextRemes

R CMD INSTALL climextRemes_${CEVER}.tar.gz

python ~/Desktop/test.py

### regular PyPI

VE=~/Desktop/ce_${CEVER}_${PYVER}
virtualenv ${VE}
source ${VE}/bin/activate

pip install climextremes==${CEVER}

### Conda from cascade channel

CHANNEL=cascade

conda create -y --name ce_${CEVER}_${PYVER} python=${PYVER}
source activate ce_${CEVER}_${PYVER}

conda install -y -c ${CHANNEL} climextremes=${CEVER}

## this will use old R climextRemes
python ~/Desktop/test.py


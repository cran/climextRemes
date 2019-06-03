#!/bin/bash

#$PYTHON setup.py install

# Add more build steps here, if they are necessary.
# $PYTHON -c "import climextremes"
# echo "install.packages(\"${pkg}\", repos=\"https://cran.rstudio.com\")" | R --no-save
#echo "install.packages(\"climextRemes\", repos=\"https://cran.r-project.org\")" | $PREFIX/bin/R --no-save
#$PREFIX/bin/pip install tzlocal

$PYTHON setup.py install

#$PYTHON -c "import climextremes"

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.

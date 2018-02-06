Python wrapping of R climextRemes package

To install climextRemes (R and Python version)

R:

install.packages("climextRemes")

To install Python version:

pip install rpy2

Depending on if you have access to standard package install directories:

python setup.py install

or

python setup.py install --prefix=/tmp/foo

export PYTHONPATH=/tmp/foo/lib/python3.5/site-packages:${PYTHONPATH}

To make documentation:

pip install sphinx

pip install sphinx-autodoc-annotation

cd docs

make html

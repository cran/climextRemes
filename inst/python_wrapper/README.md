Python wrapping of R climextRemes package
------------------------------------------

## Overview on how the Python package works

  1. The climextremes package is just the `python_wrapper/climextremes/__init__.py` file. This gets run when package is loaded.
  2. Next is an attempt to install the climextRemes R package, which should pull in the R package dependencies. Note that we try to get the version of climextRemes (R) that matches version of climextremes (python), but there is a fall-back.
  3. assuming climextRemes installs, Python calls out to R to figure out where climextRemes_wrapper is and then imports everything from climextRemes_wrapper
  4. climextRemes_wrapper goes through and creates Python fxns to wrap the R functions.

## Installation

See user-facing installation instructions in top-level README file of climextremes repository.

## Documentation

Probably best to do this in a Conda environment.

To make documentation:

  1. Make sure climextremes (Python package, presumably using setup.py method) is installed first.
  2. Put new version number into `python_wrapper/docs/build_sphinx_docs.sh` (ideally I'd figure out how to leverage `version.py`).
  3. One needs the Ubuntu `latexmk` package which is available for 18.04 but not for 16.04.
  4. If Roxygen documentation in R package has been updated, make sure you've locally installed the up-to-date R package so don't get old help info from CRAN version of R package.

```
pip install sphinx
pip install sphinx-autodoc-annotation
cd python_wrapper/docs
python make_docs.py
```

The `make_docs.py` script should create a documentation template and call `build_sphinx_docs.sh`. This script should then create the makefile and call make html and make latexpdf.

Note there are a couple warnings about formatting in fit_gev and fit_pot that I haven't been able to figure out.

```
/accounts/gen/vis/paciorek/.local/lib/python3.7/site-packages/climextremes/__init__.py:docstring of climextremes.fit_gev:84: WARNING: Block quote ends without a blank line; unexpected unindent.
/accounts/gen/vis/paciorek/.local/lib/python3.7/site-packages/climextremes/__init__.py:docstring of climextremes.fit_pot:100: WARNING: Block quote ends without a blank line; unexpected unindent.
```

## Building for installation via pip and conda.

See the README.md files in the pip and conda directories.

## Other notes

To update version number one should be able to change this in only python_wrapper/version.py.

The Python package should be OS- and Python3-version- independent. See more details here: https://www.anaconda.com/condas-new-noarch-packages/

To build the Python package see inst/pip/README.md or inst/conda/README.md.

When debugging, if needed one can directly modify contents of .py files in the installed R climextRemes package (e.g., ~/R/x86.../3.5/climextRemes/*)

numpy,pandas,rpy2,tzlocal apparently required in conda meta.yaml build: stanza. Not clear why but if try to only have python and setuptools in build: you get this:

```
Processing dependencies for climextremes==0.2.1rc8
Searching for tzlocal
Traceback (most recent call last):
  File "setup.py", line 75, in <module>
    install_requires=['numpy', 'pandas', 'rpy2', 'tzlocal'],
  File "/accounts/gen/vis/paciorek/.conda/envs/climextremes_build_min/conda-bld/climextremes_1556215226148/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold/lib/python3.7/site-packages/setuptools/__init__.py", line 145, in setup
    return distutils.core.setup(**attrs)
```

Don't need compilers in dependencies as conda r package pulls them in, specific for an OS.

noarch version seems to be Python version independent (at least for Python 3).

conda R packages always installed from source.

r-rcpparmadillo avoids an issue with installing RcppArmadillo at run-time. RcppArmadillo is pulled in by climextRemes dependencies.

rpy2 pulls in Rcpp, which is needed for a climextRemes dependency.

Could try to set up build.sh and bld.bat to install R packages (including climextRemes) at install time, but might be harder to debug user issues.

When import climextremes, the code in `__init__.py` tries to install climextRemes. It can fail (e.g. I've seen download from cran.r-project.us fail) and in that case in line 91 of __init__.py it can be unable to find climextRemes_wrapper because `__climextRemes_python_path__` not set so it can't find it in the python dir in the installed R package.

Note that the code in climextRemes_wrapper.py could be in the Python package but as it is now the same Python wrapper can point to different R versions and we bundle the Python code with the R version.

The only thing permanent in inst/python_wrapper is `climextremes/__init__.py` (the actual Python package code), docs, setup.py and metadata files (COPYRIGHTS.txt, LICENSE.txt, version.py, README.md) . All the rest has been created during experiments with building the package.

Building package uses setup.py, which is the build script for setuptools.

conda-forge approach should be able to install from pypi tarball via staged recipe.
conda-forge recipe should be able to avoid usage of build.sh/bld.bat

In some cases if one doesn't force rpy2 >=2.9.4, Conda installation will pull in rpy2=2.8.5, r=3.3.1, and there will be problems. Not clear why. 

As of August 2020, default Conda channel has rpy2 2.9.4 but pulls in pandas >= 1.0.0, which causes a run-time error because of an apparent rpy2-pandas incompatibility. conda-forge has new rpy2 but it doesn't seem possible to force use of a channel in the conda meta.yaml, nor does asking for rpy2 >=3.0.0 work when building the Conda package.

As of July 2022, I've modified the Conda deploy.sh to use `-c conda-forge` which allows me to build based on newer rpy2 and pandas, so removing the version specifications in meta.yaml.

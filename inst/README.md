See python_wrapper/README.md for information on the climextremes Python package that wraps the climextRemes R package.

See {pip,conda}/README.md for information on creating the Pip and Conda packages.

Workflow for updating climextRemes as a whole:

- Update package source and test R package via winbuilder/rhub.
- Run R tests via `test_package('climextRemes')`.
- Update setup.py to include newer versions of Python.
- Note that you need to have the local R climextRemes be the updated package with the new number while Python version will include 'rcX'. (This can cause the error message 'Invalid comparison operator in dependency: p(0.2.3rc2)' when Python climextremes tries to install R climextRemes when it is imported.)
- Rebuild Python docs via sphinx if necessary.
- Put pip version in test PyPi as rc1 and test (see test_installation.sh).
- Put Conda version in paciorek channel as rc1 and test.
- Put rc1 versions in (non-test) PyPI and conda cascade channel and test dependencies (probably not necessary if not modifying anything that doesn't affect dependencies).
- Update on CRAN making sure Python version is correct (not tagged with 'rc').
- Assuming all is well, bump to final version on (non-test) PyPI and Conda cascade channel.
- Create release on bitbucket.
- Tag release on bitbucket and update master branch.
- Post new version to Zenodo.

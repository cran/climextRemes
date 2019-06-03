Instructions for getting climextremes on PyPI.

See https://packaging.python.org/tutorials/packaging-projects.

TestPyPI for climextremes is https://test.pypi.org/manage/project/climextremes/
Note that every time you change the package seem to need to create new
version number for TestPyPI to accept the upload.

Also you need to have the twine package installed.

## Build package

```
cd python_wrapper
python setup.py sdist bdist_wheel
```

That should generate whl and .tar.gz files in `dist`.

Note that this creates an egg file but not sure why. The upload to PyPI does not seem to upload the egg.

## Upload to TestPyPI

```
# python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
python -m twine upload -u paciorek -p PASSWORD --repository-url https://test.pypi.org/legacy/ dist/*
```

If don't provide password on command line, you get "RuntimeError: No recommended backend was available. Install the keyrings.alt package if you want to use the non-recommended backends. See README.rst for details."

Note that it can be a hassle to upload new copies of files of the same version if want to test out changes. Try renaming the version number and make sure you do dist/filename rather than dist/*


## Upload to PyPI

Not sure if include '/legacy' or even if need --repository-url

```
python -m twine upload -u paciorek -p PASSWORD  dist/*
```

## Testing installation

```
# test PyPI
pip install --user --index-url https://test.pypi.org/simple/ climextremes
# regular PyPI
pip install --user climextremes
```

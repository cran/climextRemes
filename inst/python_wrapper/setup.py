"""
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
# with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
long_description = """
Functions for fitting GEV and POT (via point process fitting)
models for extremes in climate data, providing return values, return
probabilities, and return periods for stationary and nonstationary models.
Also provides differences in return values and differences in log return
probabilities for contrasts of covariate values. Functions for estimating risk
ratios for event attribution analyses, including uncertainty. Under the hood,
many of the functions use functions from 'extRemes', including for fitting the
statistical models.
"""

version = {}
with open('version.py') as f:
    exec(f.read(), version)
    
    
setup(
    name='climextremes',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=version['__version__'], 

    description='climextremes wrapper module',
    long_description=long_description,

    # The project's main homepage.
    url='https://bitbucket.org/lbl-cascade',

    # Author details
    author='Christopher Paciorek',
    author_email='paciorek@berkeley.edu',

    # Choose your license
    license='BSD License',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        #'Development Status :: 3 - Alpha',
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        #'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        #'License :: OSI Approved :: MIT License',
        'License :: OSI Approved :: BSD License',
        
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
	'Programming Language :: Python :: 3.8',
	'Programming Language :: Python :: 3.9',
	'Programming Language :: Python :: 3.10',
    ],

    # What does your project relate to?
    keywords=['statistics', 'GEV', 'POT', 'detection and attribution'],
    packages=['climextremes'],
    install_requires=['numpy', 'pandas', 'rpy2', 'tzlocal'],
)

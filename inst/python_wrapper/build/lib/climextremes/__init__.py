name = "climextremes"

import pkg_resources  # part of setuptools
__version__ = pkg_resources.require("climextremes")[0].version

def reinstall_climextremes(force, version):
    import os
    import sys
    import rpy2.robjects
    import warnings
    from rpy2.rinterface import RRuntimeWarning
    warnings.filterwarnings("ignore", category=RRuntimeWarning)

    def install_cran_default(repos = 'https://cran.us.r-project.org'):
        # sometimes failure occurs with URL issue with cran.r-project.org, so try a mirror
        try:     # this may fail with error or fail but only issue a warning
            rpy2.robjects.r("""install.packages('climextRemes',repos='{0}')""".format(repos))
            rpy2.robjects.r('''library(climextRemes)''')
        except:  
            rpy2.robjects.r("""install.packages('climextRemes',repos='https://cran.cnr.berkeley.edu')""")
            rpy2.robjects.r('''library(climextRemes)''')

    def get_devtools(repos = 'https://cran.us.r-project.org'):
        try:
            rpy2.robjects.r('''library(devtools)''')
            return True
        except:
            try:
                rpy2.robjects.r("""install.packages('devtools',repos='{0}')""".format(repos))
            except:
                rpy2.robjects.r("""install.packages('devtools',repos='https://cran.cnr.berkeley.edu')""")
        try:
            rpy2.robjects.r('''library(devtools)''')
        except:
            print("Unable to install R package devtools; needed for version-specific installation of R climextRemes.")
            return False
        return True
                

    def install_cran_specific(version, repos = 'https://cran.us.r-project.org'):
        check = get_devtools()
        # When using Conda with R provided by Conda, TAR is incorrectly set to /bin/gtar, which causes errors.
        # Simply unsetting seems to avoid the problem and not prevent installation of R packages.
        rpy2.robjects.r("""Sys.unsetenv('TAR')""")
        if check:
            try:
                rpy2.robjects.r("""install_version('climextRemes','{0}',repos='{1}', upgrade = 'never')""".format(version, repos))
                rpy2.robjects.r('''library(climextRemes)''')
            except:
                rpy2.robjects.r("""install_version('climextRemes','{0}',repos='https://cran.cnr.berkeley.edu', upgrade = 'never')""".format(version, repos))
                rpy2.robjects.r('''library(climextRemes)''')
    
    if force:
        try:
            if version is None:
                install_cran_default()
            else:
                try:
                    install_cran_specific(version)
                except:
                    print("Installation of R climextRemes version: " + version + " failed, likely because the version is not available on the CRAN R package archive.")
                    return False
        except:
            return False
    else:
        try:
            test_import = rpy2.robjects.r('''library(climextRemes)''')
            iv = rpy2.robjects.r('''packageVersion("climextRemes")''')[0]
            installed_version = str(iv[0]) + "." + str(iv[1]) + "." + str(iv[2])

            if version is not None and version != installed_version:
                print("Current version: "  + installed_version + " does not match requested version: " + version + ". Attempting installation of R climextRemes package (this may take a few minutes) ...")
                try:
                    install_cran_specific(version)
                except:
                    print("Installation of R climextRemes version: " + version + " failed (likely because the version is not on the CRAN R package archive).\nFalling back to installed version.")
        except:
            # initial import failed
            print("Attempting installation of R climextRemes package and its dependencies (this may take a few minutes) ...")
            try:
                if version is None:
                    install_cran_default()
                else:
                    try:
                        install_cran_specific(version)
                    except:
                        print("Installation of version: " + version + " failed (likely because the version is not on the CRAN R package archive).\nFalling back to installing default CRAN package.")
                        install_cran_default()
            except:
                return False
    return True

def __wrap_import():
    import os
    import sys
    import rpy2.robjects

    global __version__

    if not reinstall_climextremes(False, __version__):
            print("Installation of climextRemes failed. Please manually install climextRemes using CRAN R package archive.")
            return

    # Force R warnings to come through to Python
    rpy2.robjects.r('options(warn=1)')
    import warnings
    from rpy2.rinterface import RRuntimeWarning
    warnings.filterwarnings("always", category=RRuntimeWarning)

    __climextRemes_home__ = rpy2.robjects.r('''
    library(climextRemes)
    sp <- searchpaths()
    climextremes_path <- sp[grep("climextRemes", sp)]
    ''')[0]
    __climextRemes_python_path__ = __climextRemes_home__ + "/python"
    sys.path.append(__climextRemes_python_path__)

def __cleanup_import():
    import sys
    del sys.path[-1]
    del sys.path[-1]

__wrap_import()
from climextRemes_wrapper import *
__cleanup_import()


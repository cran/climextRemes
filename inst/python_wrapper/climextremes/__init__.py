import pkg_resources  # part of setuptools
__current_version__ = pkg_resources.require("climextremes")[0].version

# __current_version__ = "0.1.3"

def reinstall_climextremes(force, version):
  import os
  import sys
  import rpy2.robjects
  import warnings
  from rpy2.rinterface import RRuntimeWarning
  warnings.filterwarnings("ignore", category=RRuntimeWarning)

  if force:
    try:
      if version is None:
        rpy2.robjects.r("""install.packages('climextRemes',repos='http://cran.us.r-project.org')""")
      else:
        try:
          rpy2.robjects.r("""library(devtools);install_version('climextRemes','{0}',repos='http://cran.us.r-project.org')""".format(version))
        except:
          print("Installation of version: " + version + " failed, trying default CRAN package")
          rpy2.robjects.r("""install.packages('climextRemes',repos='http://cran.us.r-project.org')""")

    except:
      return False
  else:
    try:
      test_import = rpy2.robjects.r('''library(climextRemes)''')
      iv = rpy2.robjects.r('''packageVersion("climextRemes")''')[0]
      installed_version = str(iv[0]) + "." + str(iv[1]) + "." + str(iv[2])

      if version is not None and version != installed_version:
        try:
          rpy2.robjects.r("""library(devtools);install_version('climextRemes', '{0}', repos='http://cran.us.r-project.org', verbose=FALSE)""".format(version))
        except:
          print("Installation of version: " + version + " failed, trying default CRAN package")
          rpy2.robjects.r("""install.packages('climextRemes',repos='http://cran.us.r-project.org')""")

    except:
      # initial import failed
      print("Attempting installation of R climextRemes (initial time takes longer...)")

      try:
        if version is None:
          rpy2.robjects.r("""install.packages('climextRemes',repos='http://cran.us.r-project.org', verbose=FALSE)""")
        else:
          try:
            rpy2.robjects.r("""library(devtools);install_version('climextRemes','{0}',repos='http://cran.us.r-project.org')""".format(version))
          except:
            print("Installation of version: " + version + " failed, trying default CRAN package")
            rpy2.robjects.r("""install.packages('climextRemes',repos='http://cran.us.r-project.org')""")

      except:
        return False
  return True

def __wrap_import():
  import os
  import sys
  import rpy2.robjects

  global __current_version__

  if not reinstall_climextremes(False, __current_version__):
      print("Installation of climextRemes failed. Please manually install climextRemes using CRAN")
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


def __wrap_import():
  import os
  import sys
  import rpy2.robjects

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


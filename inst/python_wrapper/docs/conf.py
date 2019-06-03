import sys
import os

#sys.path.insert(0, os.path.abspath('../src/'))

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary']
source_suffix = '.rst'
master_doc = 'index'
project = u'climextremes Documentation'
copyright = u'The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) and the University of California, Berkeley'
exclude_patterns = ['_build']
pygments_style = 'sphinx'
html_theme = 'default'
autoclass_content = "both"
autodoc_default_flags = ['members']
autosummary_generate = True

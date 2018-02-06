import sys
import os

#sys.path.insert(0, os.path.abspath('../src/'))

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary']
source_suffix = '.rst'
master_doc = 'index'
project = u'climextRemes Documentation'
copyright = u'Chris Paciorek'
exclude_patterns = ['_build']
pygments_style = 'sphinx'
html_theme = 'default'
autoclass_content = "both"
autodoc_default_flags = ['members']
autosummary_generate = True

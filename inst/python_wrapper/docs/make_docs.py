import os
import types
import climextremes
#sphinx, sphinx-autodoc-annotation

x = [climextremes.__dict__.get(a) for a in dir(climextremes)
  if isinstance(climextremes.__dict__.get(a), types.FunctionType)]

output = "climextRemes Documentation\n"
output += "===========================\n\n"

output += ".. toctree::\n"
output += "  :maxdepth: 2\n"
output += "  :caption: Contents:\n\n"

#.. automodule:: climextRemes_wrapper
#  :members:

# .. autofunction:: climextremes.fit_gev

for a in x:
  if not a.__name__.startswith("__"):
    output += ".. autofunction:: climextremes.%s " % a.__name__ + "\n"
    # output += ".. autodoxymethod:: climextremes.%s " % a.__name__ + "\n"

output += "\nIndices and tables\n"
output += "==================\n\n"

output += "* :ref:`genindex`\n"
output += "* :ref:`modindex`\n"
output += "* :ref:`search`\n"

with open("index.rst", "w") as wf:
  wf.write(output)

os.system("/bin/bash build_sphinx_docs.sh")

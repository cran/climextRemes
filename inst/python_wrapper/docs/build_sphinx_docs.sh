MYDOCPATH=$PWD

DOCDIR=/tmp/climextRemes_sphinx_docs

rm -fR $DOCDIR
mkdir -p $DOCDIR

cd $DOCDIR

sphinx-quickstart -p "climextRemes Documentation" -a "Chris Paciorek" --ext-autodoc --makefile --suffix=".rst" --master="index" --sep --dot=_ -v 1.2.0 -r 1.2.0 -l en --epub --ext-doctest --ext-intersphinx --ext-todo --ext-coverage --ext-imgmath --ext-mathjax --ext-viewcode --ext-githubpages --ext-ifconfig --makefile -m --no-batchfile .

cp $MYDOCPATH/index.rst ./source/index.rst

make html
make latexpdf

cp -R ./build/latex/climextRemesDocumentation.pdf $MYDOCPATH

cp -R ./build/html $MYDOCPATH

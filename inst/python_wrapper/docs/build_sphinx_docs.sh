MYDOCPATH=$PWD

DOCDIR=/tmp/climextRemes_sphinx_docs

rm -fR $DOCDIR
mkdir -p $DOCDIR

cd $DOCDIR

sphinx-quickstart -p "climextremes Documentation" -a "Christopher Paciorek" --ext-autodoc --makefile --suffix=".rst" --master="index" --sep --dot=_ -v 0.3 -r 0.3.1 -l en --epub --ext-doctest --ext-intersphinx --ext-todo --ext-coverage --ext-imgmath --ext-mathjax --ext-viewcode --ext-githubpages --ext-ifconfig --makefile -m --no-batchfile .

cp $MYDOCPATH/index.rst ./source/index.rst

make html
make latexpdf

cp -R ./build/latex/climextremes.pdf $MYDOCPATH/climextremesDocumentation.pdf

cp -R ./build/html $MYDOCPATH

0. Ideally create a Conda environment for this; this is done currently in deploy.sh, which also shows what dependencies there are for doing the building.

1. Run code in bash deploy.sh; this should put the built package in /tmp/cebuild/noarch

2. Authenticate to Anaconda:

For use of test 'paciorek' channel
anaconda login --username paciorek --password PASSWORD

For use of real 'cascade' channel: provide password at prompt as it has a character that messes with bash passing it from the commandline
anaconda login --username cascade

3. anaconda upload /tmp/cebuild/noarch/climextremes-0.2.1-pyr351_0.tar.bz2 --force (the name should be output from 'conda build' line in deploy.sh)

4. to test:
conda install -c paciorek climextremes
conda install -c cascade climextremes


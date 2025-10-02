# ELL1Rmodel
A model version of the ELL1 model in TEMPO2, in which the Romer delay is calculated without low-eccentricity approximation

# usage
If the EDOT parameter is present in the par file, the Romer delay is calculated without the low-eccentricity approximation; otherwise, the original model is used.

# install
Replace origin ELL1model.C, then run "make && make install"

# model description
A modified version of the ELL1 model in TEMPO2 [1], in which the Romer delay is calculated without low-eccentricity approximation

# usage
If the EDOT parameter is present in the par file, the Romer delay is calculated without the low-eccentricity approximation; otherwise, the original model is used. BINARY paramter in the par file is still "ELL1".

# install
Replace original ELL1model.C, then run "make && make install"

[1] Hobbs G.~B., Edwards R.~T., Manchester R.~N., 2006, MNRAS, 369, 655. Available in https://sourceforge.net/projects/tempo2/

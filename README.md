# LGCA-persistent-nematic
C++ code for simulating an LGCA model of nematically-aligning, persistent particles, used for the paper "Individual particle persistence antagonizes global ordering in populations of nematically aligning self-propelled particles"

To compile the code, please be sure to change the path to the random number generator (RNG) and the graphical output tool (CIMG) appropriately.

Parameter values for lattice size need to be changed directly within the code before compiling; other parameter values need to be input alongside the run command in the appropriate order, e.g. ./output 500 1000 0.2 0.1 0.1 (see code atoi and atof commands) for the order of parameter input).

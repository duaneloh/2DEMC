Package function:
	Simulates data and reconstructions of a simplified version of cryptotomography. 

Included in this package:
	make_data.c
	EMC.c
	runEMC.py
	contrast.dat

Dependencies of this package:
	to compile the C files: 
		gcc, standard C libraries (stdio, stdlib, time, math)  
	to run Python file: 
		python2.7 with matplotlib, numpy, scipy (h5py optional) installed.
		For Windows users, you might want to download Enthought's python distribution (free for academic use).

How to run this package:
	1. Compile the C files into binaries:
		./runEMC.py -c
	2. Pythonic+C  mode (recommended):
		run python example 1 or 2:
		./runEMC.py -e 1 or ./runEMC.py -e 2
	3. C mode (Consult runEMC.py on how to customize runs):
		./make_data 
		./EMC

How to cite work related to this package:
	TODO: N.D. Loh, Phil. Trans. B. ....

This package was written by N. Duane Loh on Tue Jun 10 2014

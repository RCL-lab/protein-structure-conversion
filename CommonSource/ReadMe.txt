This folder contains an auxiliray file to erite the output of internal to cartesian into PDB like file, which maps our program indexing to the pdb index.
	-ItoC.h : This is the header file you should include it in the internal to cartesian CPU or GPU code.
	-writeToPDB.cpp : This is the implementation of function to write the output into pdb format.
		          This step is already existed in the CPU and GPU code but commented you can uncommented it to call the function, 
			   then get the easily to read pdb format.

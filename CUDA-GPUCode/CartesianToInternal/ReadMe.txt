The folder CartesianToInternal under CPU-GPUCode:
provides the parallel implementation of forward conversion on GPU.
parallelization of computation of internal coordinates from Cartesian coordinates, is easily paralliziable.

	- CartesiantoInternal.cu:
	   The module read the Internal coordinate of a protein and the list of atoms that make bond, angle and dihedral together. 
	   This converts Cartesian coordinates parallely to Internal coordinate by using formula. 
	   The input files of this program is trajectory.txt, bond.txt, angle.txt, improper.txt, proper.txt. 
	   The output is 4 different files for each internal Cartesian coordinate.
	   The # of atoms and # time frame can be changed based on the input. for the lysozyme input, #of atoms and time frames
	   can be kept unchanged to the value already set in the program.
				
								   

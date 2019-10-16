The folder CartesianToInternal under CPUCode-serial:
provides serial computation of internal coordinates by using the Cartesian coordinates, 
this version is only developed for the purpose of comparing between serial and our parallel method.

	- CartesiantoInternal_CPU.cpp:
	   The module read the Cartesian coordinate of a protein and the list of atoms that make bond, angle and dihedral together. 
	   Serially convert it to Internal coordinate by using formula. The input files of this program is trajectory.txt, bond.txt, 	 	    angle.txt, improper.txt, proper.txt. The output is 4 different output files for each internal Cartesian coordinate.
	   The # of atoms and # time frame can be changed based on the input with the lysozyme input we can keep the value that is already  	       set in the program.
				
								   

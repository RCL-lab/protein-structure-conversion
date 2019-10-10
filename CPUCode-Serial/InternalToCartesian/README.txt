The folder InternaltoCartesian under CPUCode-serial:
provides simple walk over the protein chain calculation using serial computation, 
this version is only developed for the purpose of comparing between serial and our parallel method.

	- InternaltoCartesian_CPU.cpp:
				The module read the intenal coordinate of a protein and serially convert it to cartesian coordinate by walking on the protein chain.
				This is the serial CPU version of the program. The input files of this program is Bond, Angle, Dihd, imprANGLE_SIDE,BOND_SIDE, ANGLE_SIDE
				DIHD_SIDE,resname. the output can be two raw data file, backboneXYZ, sideXYZ. 
				by including ItoC.h header file, we can call writeToPDB.cpp and get the pdb like output. (you only need to uncomment that part)
				The output is simplepdb.txt.
								   
	- ItoC.h:  This the header file for Internal to Cartesian coordinate, 
			   and helps write the output of internal to Cartesian conversion into pdb similar format.
               include this header internal to cartesian Source program.

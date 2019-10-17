The folder InternaltoCartesian under CUDA-GPUCode:
provides our novel parallel methods on GPU for reverse conversion on protein chain.  


	- InternaltoCartesian.cu:
	  The module read the intenal coordinate of a protein and then by breaking the residue parallely convert it to 
	  Cartesian coordinate. 
	  The input files of this program is Bond, Angle, Dihd,imprANGLE_SIDE,BOND_SIDE, ANGLE_SIDE, DIHD_SIDE,resname. 
	  The output can be two raw data file, backboneXYZ, sideXYZ. 
	  by including ItoC.h header file, we can call writeToPDB.cpp and get the pdb like output. (you only need to uncomment that part)
	  The output is simplepdb.txt.
								   
	- ItoC.h:  This the header file for Internal to Cartesian coordinate, and helps write the output of internal to Cartesian 			   conversion into pdb similar format. one can Include this header to InternaltoCartesian Source code.

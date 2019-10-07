PDB Example : 
    This folder contains of different protein pdb and some psf examples. That can be used as an input (after doing the VMD preprocessing) to the protein conversion program.
	- 1I10.pdb, 1I10.psf : https://www.rcsb.org/structure/1I10 
	  HUMAN MUSCLE L-LACTATE DEHYDROGENASE M CHAIN, 
	  These .pdb and .psf includes 8 chains, it needs to first split the chains and then extract the cordinates.	  
	- 1I10-splitchain : all 8 chains splitted from 1I10.pdb, This folder contains all the 8 chains. 
	
	
	
	- 1TIT.pdb, 1TIT.pdb: http://www.rcsb.org/structure/1TIT
   	  TITIN, This is only single chain protein of TITIN
	  
	  
	- 5rsa.pdb, 5rsa.psf:  https://www.rcsb.org/structure/5RSA
	  RIBONUCLEASE-A, This is only single chain protein, with Di-solfide bond. 
	- 5rsanowater.pdb, 5rsanowater.psf:  https://www.rcsb.org/structure/5RSA
	  RIBONUCLEASE-A with extracted water.
	  
	  
    - 6PTI.pdb : https://files.rcsb.org/view/6PTI.pdb
	  STRUCTURE OF FORM III CRYSTALS OF BOVINE PANCREATIC TRYPSIN, single chain protein.
	  
	  
	 
	-A#Number_alpha.pdb: we made these following pdb files by concatinating ALA residues to one another artificially to check the scaling of our algorithm.
		A5_alpha.pdb   : consisted of 5 ALA residues
		A10_alpha.pdb  : consisted of 10 ALA residues
		A50_alpha.pdb  : consisted of 50 ALA residues
		A100_alpha.pdb : consisted of 100 ALA residues
		A500_alpha.pdb : consisted of 500 ALA residues
		A1000_alpha.pdb: consisted of 1000 ALA residues
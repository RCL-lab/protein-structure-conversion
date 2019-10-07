VMD-Preprocess : 
	All the preprocesses with VMD Software that are needed during the protein conversion are included in this folder.
	All These files are runnable in VMD tk console, each file are provided with header and comments on how to run them. 
	Generally you need to run the them in console: Source file.tcl
	
	
	(I)   Extract-Backbone-InternalCoordinate.tcl:
	      The module extract the backbone internal coorsinates.
		  
        (II)  Extract-Sidechain-internalCordinate.tcl:
              The module extract the side chain internal coorsinates. 
    
	The output files of (I), (II) are input files for Internal to Cartesian coordinateconversion.

	
	(III)  Extract-Cartesian-Backbone&SideChain.tcl 
	      The module extract the Cartesian Coordinate of backbone and side chain of a protein. This module output only used for verification of our computation. 

	All three (I), (II) and (III) first need to load the related protein.pdb file on VMD  then Source file.tcl

	
        (IV) pdb2psf.tcl
	     This module is needed to conver a protein pdb file to psf. The psf(protein structure file) is needed when we convert the      		Cartesian coordinate to Internal coordinate, then we need to know which atoms struct bonds and angles. 
		  
		  
	(V) SplitChain.tcl
	   This module might not be needed for all proteins, but if a protein consists of multiple chain, first we need to break it into 	    multiple pdb files.
		 
		 
        (VI) OneChainMultipleFilesPDBtoPSF.tcl
	     This module is needed when we split the chains of a protein and one chain gets fragmented pdb files (generates multiple pdb   		files for one chain) In this file we convert it into one pdb and psf file.  

Protein structure conversion includes novel parallel implementation of protein internal to Cartesian (reverse) 
coordinate conversion and vice versa (forward) on GPU. 
For modeling proteins in different conformational states, e.g. by molecular dynamics or MonteCarlo methods,
two methods of atom representation are frequently used: internal coordinates and Cartesian coordinates.

Cartesian coordinates is specified with common (x,y,z) system, and Internal coordinate describes the atom position 
using the atom's bond distances and angles. 

The conversion algorithm between these two coordinate system is applicable to protein engineering and fitting protein 
structures to experimental data.

Traditionally, converting the structures of polymers (e.g.  proteins) from internal to Cartesian coordinates has been 
performedserially,  due  to  an  inherent  linear  dependency  along  the  polymer  chain.
We  show this dependency can be removed, and we observe an order of magnitude speedup using parallel processing 
on a GPU compared to serial execution on a CPU.

This repository includes the traditional CPU version of this conversion and the novel and fast parallel GPU version, 
here is a short discription of each folder in the repository:

1. PDB_Example:	Some PDB(Protein Database Bank) files  that can be used to test our program.

2. VMD-preprocess: In our conversion implementation, we generate our own inputs for both conversions. 
				   VMD-preprocessing needs to be done on PDB file or PSF (protein structure file), 
				   to get the input format ready. This preprocessing can be used with PDB_Example or 
				   other protein files.
 
3. InputInternalToCartesian : two sets of input, ready to use for Internal to Cartesian conversion.
							  These sets of input do not need to do the preprocessing any more. 

4. InputCartesianToInternal:  one set of input, ready to use for Cartesian to internal conversion.
							  This set of input do not need to do the preprocessing any more. 

5. CUDA-GPUCode: The source code for our novel parallel method. this includes both forward and reverse implementation.

6. CPUCode-Serial: The traditional serial method walking over protein structure. 
				   This implemetation is provided to compare our method with the traditional algorithm. 
				   The CPU code also includes both forward and reverse implementation.

7. CommonSource:  	The CommonSource, is needed both for CPU and GPU implementation of internal to Cartesian conversion.
					This function can be called and post-process the output and put it in the similar conventional 
					PDB format, that is easy to read.

6. VMD-CartesianCoordinateOutput: Cartesian coordinate is measured and provided to verify the accuracy output of 
								  in internal to Cartesian conversion.
								  The Cartesian output is provided for both sets of input placed in InputInternalToCartesian.






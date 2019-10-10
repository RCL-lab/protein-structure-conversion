#ifndef _ITOC_H_
#define _ITOC_H_

#include <stdlib.h>
#include <stdio.h>                         
#include <iostream>
#include <cstdlib>                                              
#include <time.h>
#include<fstream>
#include <string>
#include <iomanip> 
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cuComplex.h"
using namespace std; 
#define pi          3.141592653589793238462643383279502884L

void protein_file(double* finalxyz,double* sidexyz);

#endif
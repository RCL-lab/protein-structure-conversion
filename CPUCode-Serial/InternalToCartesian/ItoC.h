/*####################################################################################################
## FileName:    ItoC.h
## Description: This the header file for Internal to Cartesian coordinate, and helps write the output of internal to Cartesian conversion into pdb similar format.
##				include this header internal to cartesian Source program
## 				
## Author: Mahsa Bayati
##
## MIT License
## Copyright (c) 2019 RCL-lab
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
####################################################################################################
*/

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
using namespace std; 
#define pi          3.141592653589793238462643383279502884L

struct bond{
	double NCA;
	double CAC;
	double CO;
	double CN;
};

struct angle{
	double NCAC;
	double CACO;
	double CACN;
	double CNCA;
};

struct dihd{
	double NCACO;
	double NCACN;
	double CACNCA;
	double CNCAC;
};



void protein_file(double* finalxyz,double* sidexyz);
#endif
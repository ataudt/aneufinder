//aneufinder - An R-package for CNV detection in whole-genome single cell sequencing data
//Copyright (C) 2015  Aaron Taudt
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.



#ifndef UTILITY_H
#define UTILITY_H

#include <exception> // error handling
#include "logging.h" // FILE_LOG() capability
#include <R.h> // Calloc() etc.
#include <algorithm> // max_element

/* custom error handling class */
// extern statement to avoid 'multiple definition' errors
extern class exception_nan: public std::exception
{
  virtual const char* what() const throw()
  {
    return "nan detected";
  }
} nan_detected; // this line creates an object of this class

/* helpers for memory management */
double** allocDoubleMatrix(int rows, int cols);
void freeDoubleMatrix(double** matrix, int rows);
double** CallocDoubleMatrix(int rows, int cols);
void FreeDoubleMatrix(double** matrix, int rows);
int** allocIntMatrix(int rows, int cols);
void freeIntMatrix(int** matrix, int rows);
int** CallocIntMatrix(int rows, int cols);
void FreeIntMatrix(int** matrix, int rows);
bool** allocBoolMatrix(int rows, int cols);
void freeBoolMatrix(bool** matrix, int rows);
bool** CallocBoolMatrix(int rows, int cols);
void FreeBoolMatrix(bool** matrix, int rows);
double*** alloc3Ddouble(int dim1, int dim2, int dim3);
void free3Ddouble(double*** array, int dim1, int dim2);
double Max(double *a, int N);
int argMax(double *a, const int N); //return an index of the maximum (the first one if tight happens)
int intMax(int *a, int N);
int argIntMax(int *a, const int N);
double MaxMatrix(double**, int N, int M);
int MaxIntMatrix(int**, int N, int M);
double MaxDoubleMatrix(double**, int N, int M);
void printDoubleAsBinary(double someDouble);

#endif // UTILITY_H

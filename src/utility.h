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

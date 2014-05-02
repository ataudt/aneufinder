#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <cstdlib> // calloc() etc.
#include <cstring> // for memcpy() in printDoubleAsBinary()

/* helpers for memory management */
double** allocDoubleMatrix(int rows, int cols);
void freeDoubleMatrix(double** matrix, int rows);
int** allocIntMatrix(int rows, int cols);
void freeIntMatrix(int** matrix, int rows);
double*** alloc3Ddouble(int dim1, int dim2, int dim3);
void free3Ddouble(double*** array, int dim1, int dim2);
bool** allocBoolMatrix(int rows, int cols);
void freeBoolMatrix(bool** matrix, int rows);
double Max(double *a, int N);
int argMax(double *a, const int N); //return an index of the maximum (the first one if tight happens)
int intMax(int *a, int N);
int argIntMax(int *a, const int N);
double MaxMatrix(double**, int N, int M);
int MaxIntMatrix(int**, int N, int M);
double MaxDoubleMatrix(double**, int N, int M);
void printDoubleAsBinary(double someDouble);

#endif // UTILITY_H

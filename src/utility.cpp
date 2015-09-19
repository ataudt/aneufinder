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



#include "utility.h"

/* helpers for memory management */
double** allocDoubleMatrix(int rows, int cols)
{
	double** matrix = (double**) calloc(rows, sizeof(double*));
	int i;
	for (i=0; i<rows; i++)
	{
		matrix[i] = (double*) calloc(cols, sizeof(double));
	}
	return(matrix);
}

void freeDoubleMatrix(double** matrix, int rows)
{
	int i;
	for (i=0; i<rows; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

double** CallocDoubleMatrix(int rows, int cols)
{
	double** matrix = (double**) Calloc(rows, double*);
	int i;
	for (i=0; i<rows; i++)
	{
		matrix[i] = (double*) Calloc(cols, double);
	}
	return(matrix);
}

void FreeDoubleMatrix(double** matrix, int rows)
{
	int i;
	for (i=0; i<rows; i++)
	{
		Free(matrix[i]);
	}
	Free(matrix);
}

int** allocIntMatrix(int rows, int cols)
{
	int** matrix = (int**) calloc(rows, sizeof(int*));
	int i;
	for (i=0; i<rows; i++)
	{
		matrix[i] = (int*) calloc(cols, sizeof(int));
	}
	return(matrix);
}

void freeIntMatrix(int** matrix, int rows)
{
	int i;
	for (i=0; i<rows; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

int** CallocIntMatrix(int rows, int cols)
{
	int** matrix = (int**) Calloc(rows, int*);
	int i;
	for (i=0; i<rows; i++)
	{
		matrix[i] = (int*) Calloc(cols, int);
	}
	return(matrix);
}

void FreeIntMatrix(int** matrix, int rows)
{
	int i;
	for (i=0; i<rows; i++)
	{
		Free(matrix[i]);
	}
	Free(matrix);
}

bool** allocBoolMatrix(int rows, int cols)
{
	bool** matrix = (bool**) calloc(rows, sizeof(bool*));
	int i;
	for (i=0; i<rows; i++)
	{
		matrix[i] = (bool*) calloc(cols, sizeof(bool));
	}
	return(matrix);
}

void freeBoolMatrix(bool** matrix, int rows)
{
	int i;
	for (i=0; i<rows; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

bool** CallocBoolMatrix(int rows, int cols)
{
	bool** matrix = (bool**) Calloc(rows, bool*);
	int i;
	for (i=0; i<rows; i++)
	{
		matrix[i] = (bool*) Calloc(cols, bool);
	}
	return(matrix);
}

void FreeBoolMatrix(bool** matrix, int rows)
{
	int i;
	for (i=0; i<rows; i++)
	{
		Free(matrix[i]);
	}
	Free(matrix);
}

double*** alloc3Ddouble(int dim1, int dim2, int dim3)
{
	int i,j;
	double *** array = (double ***)malloc(dim1*sizeof(double**));

	for (i = 0; i< dim1; i++)
	{
		array[i] = (double **) malloc(dim2*sizeof(double *));
		for (j = 0; j < dim2; j++)
		{
			array[i][j] = (double *)malloc(dim3*sizeof(double));
		}
	}
	return array;
}

void free3Ddouble(double*** array, int dim1, int dim2)
{
	for (int i=0; i < dim1; i++)
	{
		freeDoubleMatrix(array[i], dim2);
	}
	free(array);
}

int argMax(double *a, const int N)
{
	double maximum=a[0];
	int argmax=0;
	for(int i=0;i<N;i++)
	{
		if(maximum<a[i])
		{
			argmax=i;
			maximum=a[i];
		}
	}
	return argmax;
}

int argIntMax(int *a, const int N)
{
	int maximum=a[0];
	int argmax=0;
	for(int i=0;i<N;i++)
	{
		if(maximum<a[i])
		{
			argmax=i;
			maximum=a[i];
		}
	}
	return argmax;
}

double Max(double *a, int N)
{
	double maximum=a[0];
	for(int i=0;i<N;i++)
	{
		if(maximum<a[i]) maximum=a[i];
	}
	return maximum;
}

int intMax(int *a, int N)
{
	int maximum=a[0];
	for(int i=0;i<N;i++)
	{
		if(maximum<a[i]) maximum=a[i];
	}
	return maximum;
}

double MaxMatrix(double **a, int N, int M)
{
	double maximum=a[0][0];
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			if(maximum<a[i][j]) maximum=a[i][j];
		}
	}
	return maximum;
}
int MaxIntMatrix(int** a, int N, int M)
{
	int maximum=a[0][0];
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			if(maximum<a[i][j]) maximum=a[i][j];
		}
	}
	return maximum;
}

double MaxDoubleMatrix(double** a, int N, int M)
{
	double maximum=a[0][0];
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			if(maximum<a[i][j]) maximum=a[i][j];
		}
	}
	return maximum;
}

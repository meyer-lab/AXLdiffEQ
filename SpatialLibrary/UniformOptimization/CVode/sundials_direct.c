/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2010/12/01 22:46:56 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for operations to be used by a
 * generic direct linear solver.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include "sundials_direct.h"
#include "sundials_nvector.h"

DlsMat NewDenseMat(long int M, long int N)
{
  DlsMat A;
  long int j;

  if ( (M <= 0) || (N <= 0) ) return(NULL);

  A = NULL;
  A = (DlsMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = (double *) malloc(M * N * sizeof(double));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }
  A->cols = (double **) malloc(N * sizeof(double *));
  if (A->cols == NULL) {
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }

  for (j=0; j < N; j++) A->cols[j] = A->data + j * M;

  A->M = M;
  A->N = N;
  A->ldim = M;
  A->ldata = M*N;

  A->type = SUNDIALS_DENSE;

  return(A);
}

double **newDenseMat(long int m, long int n)
{
  long int j;
  double **a;

  if ( (n <= 0) || (m <= 0) ) return(NULL);

  a = NULL;
  a = (double **) malloc(n * sizeof(double *));
  if (a == NULL) return(NULL);

  a[0] = NULL;
  a[0] = (double *) malloc(m * n * sizeof(double));
  if (a[0] == NULL) {
    free(a); a = NULL;
    return(NULL);
  }

  for (j=1; j < n; j++) a[j] = a[0] + j * m;

  return(a);
}


DlsMat NewBandMat(long int N, long int mu, long int ml, long int smu)
{
  DlsMat A;
  long int j, colSize;

  if (N <= 0) return(NULL);
  
  A = NULL;
  A = (DlsMat) malloc(sizeof *A);
  if (A == NULL) return (NULL);

  colSize = smu + ml + 1;
  A->data = NULL;
  A->data = (double *) malloc(N * colSize * sizeof(double));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }

  A->cols = NULL;
  A->cols = (double **) malloc(N * sizeof(double *));
  if (A->cols == NULL) {
    free(A->data);
    free(A); A = NULL;
    return(NULL);
  }

  for (j=0; j < N; j++) A->cols[j] = A->data + j * colSize;

  A->M = N;
  A->N = N;
  A->mu = mu;
  A->ml = ml;
  A->s_mu = smu;
  A->ldim =  colSize;
  A->ldata = N * colSize;

  A->type = SUNDIALS_BAND;

  return(A);
}

double **newBandMat(long int n, long int smu, long int ml)
{
  double **a;
  long int j, colSize;

  if (n <= 0) return(NULL);

  a = NULL;
  a = (double **) malloc(n * sizeof(double *));
  if (a == NULL) return(NULL);

  colSize = smu + ml + 1;
  a[0] = NULL;
  a[0] = (double *) malloc(n * colSize * sizeof(double));
  if (a[0] == NULL) {
    free(a); a = NULL;
    return(NULL);
  }

  for (j=1; j < n; j++) a[j] = a[0] + j * colSize;

  return(a);
}

void DestroyMat(DlsMat A)
{
  free(A->data);  A->data = NULL;
  free(A->cols);
  free(A); A = NULL;
}

void destroyMat(double **a)
{
  free(a[0]); a[0] = NULL;
  free(a); a = NULL;
}

int *NewIntArray(int N)
{
  int *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (int *) malloc(N * sizeof(int));

  return(vec);
}

int *newIntArray(int n)
{
  int *v;

  if (n <= 0) return(NULL);

  v = NULL;
  v = (int *) malloc(n * sizeof(int));

  return(v);
}

long int *NewLintArray(long int N)
{
  long int *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (long int *) malloc(N * sizeof(long int));

  return(vec);
}

long int *newLintArray(long int n)
{
  long int *v;

  if (n <= 0) return(NULL);

  v = NULL;
  v = (long int *) malloc(n * sizeof(long int));

  return(v);
}

double *NewRealArray(long int N)
{
  double *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (double *) malloc(N * sizeof(double));

  return(vec);
}

double *newRealArray(long int m)
{
  double *v;

  if (m <= 0) return(NULL);

  v = NULL;
  v = (double *) malloc(m * sizeof(double));

  return(v);
}

void DestroyArray(void *V)
{ 
  free(V); 
  V = NULL;
}

void destroyArray(void *v)
{
  free(v); 
  v = NULL;
}


void AddIdentity(DlsMat A)
{
  long int i;

  switch (A->type) {

  case SUNDIALS_DENSE:
    for (i=0; i<A->N; i++) A->cols[i][i] += 1.0;
    break;

  case SUNDIALS_BAND:
    for (i=0; i<A->M; i++) A->cols[i][A->s_mu] += 1.0;
    break;

  }

}


void SetToZero(DlsMat A)
{
  long int i, j, colSize;
  double *col_j;

  switch (A->type) {

  case SUNDIALS_DENSE:
    
    for (j=0; j<A->N; j++) {
      col_j = A->cols[j];
      for (i=0; i<A->M; i++)
        col_j[i] = 0.0;
    }

    break;

  case SUNDIALS_BAND:

    colSize = A->mu + A->ml + 1;
    for (j=0; j<A->M; j++) {
      col_j = A->cols[j] + A->s_mu - A->mu;
      for (i=0; i<colSize; i++)
        col_j[i] = 0.0;
    }

    break;

  }

}


void PrintMat(DlsMat A)
{
  long int i, j, start, finish;
  double **a;

  switch (A->type) {

  case SUNDIALS_DENSE:

    printf("\n");
    for (i=0; i < A->M; i++) {
      for (j=0; j < A->N; j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%12Lg  ", DENSE_ELEM(A,i,j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%12lg  ", DENSE_ELEM(A,i,j));
#else
        printf("%12g  ", DENSE_ELEM(A,i,j));
#endif
      }
      printf("\n");
    }
    printf("\n");
    
    break;

  case SUNDIALS_BAND:

    a = A->cols;
    printf("\n");
    for (i=0; i < A->N; i++) {
      start = MAX(0,i-A->ml);
      finish = MIN(A->N-1,i+A->mu);
      for (j=0; j < start; j++) printf("%12s  ","");
      for (j=start; j <= finish; j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%12Lg  ", a[j][i-j+A->s_mu]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%12lg  ", a[j][i-j+A->s_mu]);
#else
        printf("%12g  ", a[j][i-j+A->s_mu]);
#endif
      }
      printf("\n");
    }
    printf("\n");
    
    break;

  }

}



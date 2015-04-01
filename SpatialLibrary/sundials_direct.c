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

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

DlsMat NewDenseMat(size_t M, size_t N) {
  DlsMat A;

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

  for (size_t j=0; j < N; j++) A->cols[j] = A->data + j * M;

  A->M = M;
  A->N = N;
  A->ldim = M;
  A->ldata = M*N;
  A->type = SUNDIALS_DENSE;

  return(A);
}

void DestroyMat(DlsMat A)
{
  free(A->data);  A->data = NULL;
  free(A->cols);
  free(A); A = NULL;
}

long int *NewLintArray(long int N)
{
  long int *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (long int *) malloc((size_t) N * sizeof(long int));

  return(vec);
}

long int *newLintArray(long int n)
{
  long int *v;

  if (n <= 0) return(NULL);

  v = NULL;
  v = (long int *) malloc((size_t) n * sizeof(long int));

  return(v);
}

void DestroyArray(void *V)
{ 
  free(V); 
  V = NULL;
}

void AddIdentity(DlsMat A)
{
  long int i;

    for (i=0; i<A->N; i++) A->cols[i][i] += ONE;


}


void SetToZero(DlsMat A)
{
  long int i, j;
  double *col_j;
    
    for (j=0; j<A->N; j++) {
      col_j = A->cols[j];
      for (i=0; i<A->M; i++)
        col_j[i] = ZERO;
    }


}


void PrintMat(DlsMat A)
{
  long int i, j;

    printf("\n");
    for (i=0; i < A->M; i++) {
      for (j=0; j < A->N; j++) {
        printf("%12lg  ", DENSE_ELEM(A,i,j));
      }
      printf("\n");
    }
    printf("\n");


}



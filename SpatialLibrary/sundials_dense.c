/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2010/12/01 22:46:56 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic package of dense
 * matrix operations.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include "sundials_dense.h"
#include "sundials_nvector.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/*
 * -----------------------------------------------------
 * Functions working on DlsMat
 * -----------------------------------------------------
 */



long int denseGETRF(double **a, size_t m, size_t n, long int *p)
{
  long int i, j, k, l;
  double *col_j, *col_k;
  double temp, mult, a_kj;

  /* k-th elimination step number */
  for (k=0; k < n; k++) {

    col_k  = a[k];

    /* find l = pivot row number */
    l=k;
    for (i=k+1; i < m; i++)
      if (fabs(col_k[i]) > fabs(col_k[l])) l=i;
    p[k] = l;

    /* check for zero pivot element */
    if (col_k[l] == ZERO) return(k+1);
    
    /* swap a(k,1:n) and a(l,1:n) if necessary */    
    if ( l!= k ) {
      for (i=0; i<n; i++) {
        temp = a[i][l];
        a[i][l] = a[i][k];
        a[i][k] = temp;
      }
    }

    /* Scale the elements below the diagonal in
     * column k by 1.0/a(k,k). After the above swap
     * a(k,k) holds the pivot element. This scaling
     * stores the pivot row multipliers a(i,k)/a(k,k)
     * in a(i,k), i=k+1, ..., m-1.                      
     */
    mult = ONE/col_k[k];
    for(i=k+1; i < m; i++) col_k[i] *= mult;

    /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., m-1 */
    /* row k is the pivot row after swapping with row l.      */
    /* The computation is done one column at a time,          */
    /* column j=k+1, ..., n-1.                                */

    for (j=k+1; j < n; j++) {

      col_j = a[j];
      a_kj = col_j[k];

      /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)  */
      /* a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k) */

      if (a_kj != ZERO) {
	for (i=k+1; i < m; i++)
	  col_j[i] -= a_kj * col_k[i];
      }
    }
  }

  /* return 0 to indicate success */

  return(0);
}

void denseGETRS(double **a, size_t n, long int *p, double *b)
{
  size_t i, pk;
  double *col_k, tmp;

  /* Permute b, based on pivot information in p */
  for (size_t k=0; k<n; k++) {
    pk = (size_t) p[k];
    if(pk != k) {
      tmp = b[k];
      b[k] = b[pk];
      b[pk] = tmp;
    }
  }

  /* Solve Ly = b, store solution y in b */
  for (size_t k=0; k<n-1; k++) {
    col_k = a[k];
    for (i=k+1; i<n; i++) b[i] -= col_k[i]*b[k];
  }

  /* Solve Ux = y, store solution x in b */
  for (size_t k = n-1; k > 0; k--) {
    col_k = a[k];
    b[k] /= col_k[k];
    for (i=0; i<k; i++) b[i] -= col_k[i]*b[k];
  }
  b[0] /= a[0][0];

}

static void denseCopy(double ** const a, double ** const b, const size_t m, const size_t n)
{
  long int i, j;
  double *a_col_j, *b_col_j;

  for (j=0; j < n; j++) {
    a_col_j = a[j];
    b_col_j = b[j];
    for (i=0; i < m; i++)
      b_col_j[i] = a_col_j[i];
  }

}

static void denseScale(double c, double **a, size_t m, size_t n)
{
  size_t i;
  double *col_j;

  for (size_t j=0; j < n; j++) {
    col_j = a[j];
    for (i=0; i < m; i++)
      col_j[i] *= c;
  }
}

void denseAddIdentity(double **a, long int n)
{
  long int i;
  
  for (i=0; i < n; i++) a[i][i] += ONE;
}

long int DenseGETRF(DlsMat A, long int *p)
{
    return(denseGETRF(A->cols, A->M, A->N, p));
}

void DenseGETRS(DlsMat A, long int *p, double *b)
{
    denseGETRS(A->cols, A->N, p, b);
}

void DenseCopy(DlsMat A, DlsMat B)
{
    denseCopy(A->cols, B->cols, A->M, A->N);
}

void DenseScale(double c, DlsMat A)
{
    denseScale(c, A->cols, A->M, A->N);
}

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
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "sundials_dense.h"
#define TWO  (2.0)

/*
 * -----------------------------------------------------
 * Functions working on DlsMat
 * -----------------------------------------------------
 */

long int DenseGETRF(DlsMat A, long int *p) {
  return(denseGETRF(A->cols, A->M, A->N, p));
}



void DenseCopy(DlsMat A, DlsMat B) {
  denseCopy(A->cols, B->cols, A->M, A->N);
}

void DenseScale(double c, DlsMat A) {
  denseScale(c, A->cols, A->M, A->N);
}

long int denseGETRF(double **a, long int m, long int n, long int *p) {
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
    if (col_k[l] == 0.0) return(k+1);
    
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
    mult = 1.0/col_k[k];
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

      if (a_kj != 0.0) {
	for (i=k+1; i < m; i++)
	  col_j[i] -= a_kj * col_k[i];
      }
    }
  }

  /* return 0 to indicate success */

  return(0);
}



void denseCopy(double **a, double **b, long int m, long int n) {
  for (size_t j = 0; j < n; j++) memcpy(b[j], a[j], m*sizeof(double));
}

void denseScale(double c, double **a, long int m, long int n) {
  for (size_t j=0; j < n; j++) {
    for (size_t i=0; i < m; i++) a[j][i] *= c;
  }
}

void denseAddIdentity(double **a, long int n) {
  for (size_t i=0; i < n; i++) a[i][i] += 1.0;
}

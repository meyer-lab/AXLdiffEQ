/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2010/12/01 22:17:18 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for a generic package of DENSE matrix
 * operations, based on the DlsMat type defined in sundials_direct.h.
 *
 * There are two sets of dense solver routines listed in
 * this file: one set uses type DlsMat defined below and the
 * other set uses the type double ** for dense matrix arguments.
 * Routines that work with the type DlsMat begin with "Dense".
 * Routines that work with double** begin with "dense". 
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_DENSE_H
#define _SUNDIALS_DENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "sundials_direct.h"

/*
 * -----------------------------------------------------------------
 * Functions: DenseGETRF and DenseGETRS
 * -----------------------------------------------------------------
 * DenseGETRF performs the LU factorization of the M by N dense
 * matrix A. This is done using standard Gaussian elimination
 * with partial (row) pivoting. Note that this applies only
 * to matrices with M >= N and full column rank.
 *
 * A successful LU factorization leaves the matrix A and the
 * pivot array p with the following information:
 *
 * (1) p[k] contains the row number of the pivot element chosen
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.
 *
 * (2) If the unique LU factorization of A is given by PA = LU,
 *     where P is a permutation matrix, L is a lower trapezoidal
 *     matrix with all 1's on the diagonal, and U is an upper
 *     triangular matrix, then the upper triangular part of A
 *     (including its diagonal) contains U and the strictly lower
 *     trapezoidal part of A contains the multipliers, I-L.
 *
 * For square matrices (M=N), L is unit lower triangular.
 *
 * DenseGETRF returns 0 if successful. Otherwise it encountered
 * a zero diagonal element during the factorization. In this case
 * it returns the column index (numbered from one) at which
 * it encountered the zero.
 *
 * DenseGETRS solves the N-dimensional system A x = b using
 * the LU factorization in A and the pivot information in p
 * computed in DenseGETRF. The solution x is returned in b. This
 * routine cannot fail if the corresponding call to DenseGETRF
 * did not fail.
 * DenseGETRS does NOT check for a square matrix!
 *
 * -----------------------------------------------------------------
 * DenseGETRF and DenseGETRS are simply wrappers around denseGETRF
 * and denseGETRS, respectively, which perform all the work by
 * directly accessing the data in the DlsMat A (i.e. the field cols)
 * -----------------------------------------------------------------
 */

long int DenseGETRF(DlsMat A, long int *p);
void DenseGETRS(DlsMat A, long int *p, double *b);

size_t denseGETRF(double **a, size_t m, size_t n, long int *p);
void denseGETRS(double **a, long int n, long int *p, double *b);



/*
 * -----------------------------------------------------------------
 * Function : DenseCopy
 * -----------------------------------------------------------------
 * DenseCopy copies the contents of the M-by-N matrix A into the
 * M-by-N matrix B.
 * 
 * DenseCopy is a wrapper around denseCopy which accesses the data
 * in the DlsMat A and B (i.e. the fields cols)
 * -----------------------------------------------------------------
 */

void DenseCopy(DlsMat A, DlsMat B);
void denseCopy(double **a, double **b, size_t m, size_t n);

/*
 * -----------------------------------------------------------------
 * Function: DenseScale
 * -----------------------------------------------------------------
 * DenseScale scales the elements of the M-by-N matrix A by the
 * constant c and stores the result back in A.
 *
 * DenseScale is a wrapper around denseScale which performs the actual
 * scaling by accessing the data in the DlsMat A (i.e. the field
 * cols).
 * -----------------------------------------------------------------
 */

void DenseScale(double c, DlsMat A);
void denseScale(double c, double **a, size_t m, size_t n);


/*
 * -----------------------------------------------------------------
 * Function: denseAddIdentity
 * -----------------------------------------------------------------
 * denseAddIdentity adds the identity matrix to the n-by-n matrix
 * stored in the double** arrays.
 * -----------------------------------------------------------------
 */

void denseAddIdentity(double **a, size_t n);

#ifdef __cplusplus
}
#endif

#endif

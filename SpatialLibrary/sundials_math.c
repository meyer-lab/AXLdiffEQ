/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006/07/05 15:32:38 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a simple C-language math
 * library.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sundials_math.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

double RPowerI(double base, int exponent)
{
  int i, expt;
  double prod;

  prod = ONE;
  expt = abs(exponent);
  for(i = 1; i <= expt; i++) prod *= base;
  if (exponent < 0) prod = ONE/prod;
  return(prod);
}

double RPowerR(double base, double exponent)
{
  if (base <= ZERO) return(ZERO);
  return(pow(base, exponent));
}

double RSqrt(double x)
{
  if (x <= ZERO) return(ZERO);
  return(sqrt(x));
}

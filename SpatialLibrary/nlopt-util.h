/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#ifndef NLOPT_UTIL_H
#define NLOPT_UTIL_H

#include <stdlib.h>
#include <math.h>

#include "nlopt.h"

/* workaround for Solaris + gcc 3.4.x bug (see configure.ac) */
#if defined(__GNUC__) && defined(REPLACEMENT_HUGE_VAL)
#  undef HUGE_VAL
#  define HUGE_VAL REPLACEMENT_HUGE_VAL
#endif



/* Define if syscall(SYS_gettid) available. */
#define HAVE_GETTID_SYSCALL 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if the isinf() function/macro is available. */
#define HAVE_ISINF 1

/* Define if the isnan() function/macro is available. */
#define HAVE_ISNAN 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `qsort_r' function. */
#define HAVE_QSORT_R 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `time' function. */
#define HAVE_TIME 1

/* Define to 1 if the system has the type `uint32_t'. */
#define HAVE_UINT32_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to the sub-directory in which libtool stores uninstalled libraries.
 */
#define LT_OBJDIR ".libs/"

/* replacement for broken HUGE_VAL macro, if needed */
/* #undef REPLACEMENT_HUGE_VAL */

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to C thread-local keyword, or to nothing if this is not supported in
 your compiler. */
#define THREADLOCAL __thread

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Version number of package */
#define VERSION "2.4.2"

/* Define if compiled including C++-based routines */
/* #undef WITH_CXX */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
 calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif








#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

int nlopt_isinf(double x);




/* pseudorandom number generation by Mersenne twister algorithm */
extern double nlopt_urand(double a, double b);
extern int nlopt_iurand(int n);
extern double nlopt_nrand(double mean, double stddev);

/* Sobol' low-discrepancy-sequence generation */
typedef struct nlopt_soboldata_s *nlopt_sobol;
extern nlopt_sobol nlopt_sobol_create(unsigned sdim);
extern void nlopt_sobol_destroy(nlopt_sobol s);
extern void nlopt_sobol_next01(nlopt_sobol s, double *x);
extern void nlopt_sobol_next(nlopt_sobol s, double *x,
			    const double *lb, const double *ub);
extern void nlopt_sobol_skip(nlopt_sobol s, unsigned n, double *x);

/* stopping criteria */
typedef struct {
     unsigned n;
     double minf_max;
     double ftol_rel;
     double ftol_abs;
     double xtol_rel;
     const double *xtol_abs;
     int nevals, maxeval;
     double maxtime, start;
     int *force_stop;
} nlopt_stopping;
extern int nlopt_stop_f(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_ftol(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_x(const nlopt_stopping *stop, 
			const double *x, const double *oldx);
extern int nlopt_stop_dx(const nlopt_stopping *stop, 
			 const double *x, const double *dx);
extern int nlopt_stop_xs(const nlopt_stopping *stop, 
			 const double *xs, const double *oldxs,
			 const double *scale_min, const double *scale_max);
extern int nlopt_stop_evals(const nlopt_stopping *stop);
extern int nlopt_stop_evalstime(const nlopt_stopping *stop);
extern int nlopt_stop_forced(const nlopt_stopping *stop);

/* for local optimizations, temporarily setting eval/time limits */
extern nlopt_result nlopt_optimize_limited(nlopt_opt opt, 
					   double *x, double *minf,
					   int maxevals, double maxtime);

/* data structure for nonlinear inequality or equality constraint
   (f <= 0 or f = 0, respectively).  tol (>= 0) is a tolerance
   that is used for stopping criteria -- the point is considered
   "feasible" for purposes of stopping if the constraint is violated
   by at most tol. */
typedef struct {
     unsigned m; /* dimensional of constraint: mf maps R^n -> R^m */
     nlopt_func f; /* one-dimensional constraint, requires m == 1 */
     nlopt_mfunc mf;
     nlopt_precond pre; /* preconditioner for f (NULL if none or if mf) */
     void *f_data;
     double *tol;
} nlopt_constraint;

extern unsigned nlopt_count_constraints(unsigned p, const nlopt_constraint *c);
extern unsigned nlopt_max_constraint_dim(unsigned p, const nlopt_constraint *c);
extern void nlopt_eval_constraint(double *result, double *grad,
				  const nlopt_constraint *c,
				  unsigned n, const double *x);

/* rescale.c: */
double *nlopt_compute_rescaling(unsigned n, const double *dx);
double *nlopt_new_rescaled(unsigned n, const double *s, const double *x);
void nlopt_rescale(unsigned n, const double *s, const double *x, double *xs);
void nlopt_unscale(unsigned n, const double *s, const double *x, double *xs);
void nlopt_reorder_bounds(unsigned n, double *lb, double *ub);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif



#ifndef NLOPT_INTERNAL_H
#define NLOPT_INTERNAL_H

#include "nlopt.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
    
    /*********************************************************************/
    
    struct nlopt_opt_s {
        nlopt_algorithm algorithm; /* the optimization algorithm (immutable) */
        unsigned n; /* the dimension of the problem (immutable) */
        
        nlopt_func f; void *f_data; /* objective function to minimize */
        nlopt_precond pre; /* optional preconditioner for f (NULL if none) */
        int maximize; /* nonzero if we are maximizing, not minimizing */
        
        double *lb, *ub; /* lower and upper bounds (length n) */
        
        unsigned m; /* number of inequality constraints */
        unsigned m_alloc; /* number of inequality constraints allocated */
        nlopt_constraint *fc; /* inequality constraints, length m_alloc */
        
        unsigned p; /* number of equality constraints */
        unsigned p_alloc; /* number of inequality constraints allocated */
        nlopt_constraint *h; /* equality constraints, length p_alloc */
        
        nlopt_munge munge_on_destroy, munge_on_copy; /* hack for wrappers */
        
        /* stopping criteria */
        double stopval; /* stop when f reaches stopval or better */
        double ftol_rel, ftol_abs; /* relative/absolute f tolerances */
        double xtol_rel, *xtol_abs; /* rel/abs x tolerances */
        int maxeval; /* max # evaluations */
        
        int force_stop; /* if nonzero, force a halt the next time we
                         try to evaluate the objective during optimization */
        /* when local optimization is used, we need a force_stop in the
         parent object to force a stop in child optimizations */
        struct nlopt_opt_s *force_stop_child;
        
        /* algorithm-specific parameters */
        nlopt_opt local_opt; /* local optimizer */
        unsigned stochastic_population; /* population size for stochastic algs */
        double *dx; /* initial step sizes (length n) for nonderivative algs */
        unsigned vector_storage; /* max subspace dimension (0 for default) */
        
        void *work; /* algorithm-specific workspace during optimization */
    };
    
    
    /*********************************************************************/
    /* global defaults set by deprecated API: */
    
    extern nlopt_algorithm nlopt_local_search_alg_deriv;
    extern nlopt_algorithm nlopt_local_search_alg_nonderiv;
    extern int nlopt_local_search_maxeval;
    extern unsigned nlopt_stochastic_population;
    
    /*********************************************************************/
    
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* NLOPT_INTERNAL_H */


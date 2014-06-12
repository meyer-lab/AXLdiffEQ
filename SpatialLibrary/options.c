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

#include "cobyla.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "nlopt-util.h"

#define NLOPT_NUM_ALGORITHMS 1

/*************************************************************************/

void NLOPT_STDCALL nlopt_destroy(nlopt_opt opt)
{
     if (opt) {
	  unsigned i;
	  if (opt->munge_on_destroy) {
	       nlopt_munge munge = opt->munge_on_destroy;
	       munge(opt->f_data);
	       for (i = 0; i < opt->m; ++i)
		    munge(opt->fc[i].f_data);
	       for (i = 0; i < opt->p; ++i)
		    munge(opt->h[i].f_data);
	  }
	  for (i = 0; i < opt->m; ++i)
	       free(opt->fc[i].tol);
	  for (i = 0; i < opt->p; ++i)
	       free(opt->h[i].tol);
	  free(opt->lb); free(opt->ub);
	  free(opt->xtol_abs);
	  free(opt->fc);
	  free(opt->h);
	  nlopt_destroy(opt->local_opt);
	  free(opt->dx);
	  free(opt->work);
	  free(opt);
     }
}

nlopt_opt NLOPT_STDCALL nlopt_create(nlopt_algorithm algorithm, unsigned n)
{
     nlopt_opt opt;

     opt = (nlopt_opt) malloc(sizeof(struct nlopt_opt_s));
     if (opt) {
	  opt->algorithm = algorithm;
	  opt->n = n;
	  opt->f = NULL; opt->f_data = NULL; opt->pre = NULL;
	  opt->maximize = 0;
	  opt->munge_on_destroy = opt->munge_on_copy = NULL;

	  opt->lb = opt->ub = NULL;
	  opt->m = opt->m_alloc = 0;
	  opt->fc = NULL;
	  opt->p = opt->p_alloc = 0;
	  opt->h = NULL;

	  opt->stopval = -HUGE_VAL;
	  opt->ftol_rel = opt->ftol_abs = 0;
	  opt->xtol_rel = 0; opt->xtol_abs = NULL;
	  opt->maxeval = 0;
	  opt->force_stop = 0;
	  opt->force_stop_child = NULL;

	  opt->local_opt = NULL;
	  opt->stochastic_population = 0;
	  opt->vector_storage = 0;
	  opt->dx = NULL;
	  opt->work = NULL;

	  if (n > 0) {
	       opt->lb = (double *) malloc(sizeof(double) * (n));
	       if (!opt->lb) goto oom;
	       opt->ub = (double *) malloc(sizeof(double) * (n));
	       if (!opt->ub) goto oom;
	       opt->xtol_abs = (double *) malloc(sizeof(double) * (n));
	       if (!opt->xtol_abs) goto oom;
	       nlopt_set_lower_bounds1(opt, -HUGE_VAL);
	       nlopt_set_upper_bounds1(opt, +HUGE_VAL);
	       nlopt_set_xtol_abs1(opt, 0.0);
	  }
     }

     return opt;

oom:
     nlopt_destroy(opt);
     return NULL;
}

/*************************************************************************/

nlopt_result NLOPT_STDCALL nlopt_set_precond_min_objective(nlopt_opt opt,
							   nlopt_func f, 
							   nlopt_precond pre,
							   void *f_data)
{
     if (opt) {
	  if (opt->munge_on_destroy) opt->munge_on_destroy(opt->f_data);
	  opt->f = f; opt->f_data = f_data; opt->pre = pre;
	  opt->maximize = 0;
	  if (nlopt_isinf(opt->stopval) && opt->stopval > 0)
	       opt->stopval = -HUGE_VAL; /* switch default from max to min */
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result NLOPT_STDCALL nlopt_set_min_objective(nlopt_opt opt,
						   nlopt_func f, void *f_data)
{
     return nlopt_set_precond_min_objective(opt, f, NULL, f_data);
}

nlopt_result NLOPT_STDCALL nlopt_set_precond_max_objective(nlopt_opt opt, 
							   nlopt_func f, 
							   nlopt_precond pre,
							   void *f_data)
{
     if (opt) {
	  if (opt->munge_on_destroy) opt->munge_on_destroy(opt->f_data);
	  opt->f = f; opt->f_data = f_data; opt->pre = pre;
	  opt->maximize = 1;
	  if (nlopt_isinf(opt->stopval) && opt->stopval < 0)
	       opt->stopval = +HUGE_VAL; /* switch default from min to max */
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result NLOPT_STDCALL nlopt_set_max_objective(nlopt_opt opt,
						   nlopt_func f, void *f_data)
{
     return nlopt_set_precond_max_objective(opt, f, NULL, f_data);
}

/*************************************************************************/

nlopt_result
NLOPT_STDCALL nlopt_set_lower_bounds(nlopt_opt opt, const double *lb)
{
     if (opt && (opt->n == 0 || lb)) {
	  memcpy(opt->lb, lb, sizeof(double) * (opt->n));
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result
NLOPT_STDCALL nlopt_set_lower_bounds1(nlopt_opt opt, double lb)
{
     if (opt) {
	  unsigned i;
	  for (i = 0; i < opt->n; ++i)
	       opt->lb[i] = lb;
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result
NLOPT_STDCALL nlopt_get_lower_bounds(const nlopt_opt opt, double *lb)
{
     if (opt && (opt->n == 0 || lb)) {
	  memcpy(lb, opt->lb, sizeof(double) * (opt->n));
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result
NLOPT_STDCALL nlopt_set_upper_bounds(nlopt_opt opt, const double *ub)
{
     if (opt && (opt->n == 0 || ub)) {
	  memcpy(opt->ub, ub, sizeof(double) * (opt->n));
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result
NLOPT_STDCALL nlopt_set_upper_bounds1(nlopt_opt opt, double ub)
{
     if (opt) {
	  unsigned i;
	  for (i = 0; i < opt->n; ++i)
	       opt->ub[i] = ub;
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result
NLOPT_STDCALL nlopt_get_upper_bounds(const nlopt_opt opt, double *ub)
{
     if (opt && (opt->n == 0 || ub)) {
	  memcpy(ub, opt->ub, sizeof(double) * (opt->n));
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

/*************************************************************************/

#define AUGLAG_ALG(a) ((a) == NLOPT_AUGLAG ||		\
	               (a) == NLOPT_AUGLAG_EQ ||        \
	               (a) == NLOPT_LN_AUGLAG ||        \
		       (a) == NLOPT_LN_AUGLAG_EQ ||     \
		       (a) == NLOPT_LD_AUGLAG ||        \
		       (a) == NLOPT_LD_AUGLAG_EQ)





/*************************************************************************/

#define SET(param, T, arg)						\
   nlopt_result NLOPT_STDCALL nlopt_set_##param(nlopt_opt opt, T arg)	\
   {									\
	if (opt) {							\
	     opt->arg = arg;						\
	     return NLOPT_SUCCESS;					\
	}								\
	return NLOPT_INVALID_ARGS;					\
   }


#define GET(param, T, arg) T NLOPT_STDCALL	\
   nlopt_get_##param(const nlopt_opt opt) {	\
        return opt->arg;			\
   }

#define GETSET(param, T, arg) GET(param, T, arg) SET(param, T, arg)

GETSET(stopval, double, stopval)

GETSET(ftol_rel, double, ftol_rel)
GETSET(ftol_abs, double, ftol_abs)
GETSET(xtol_rel, double, xtol_rel)

nlopt_result
NLOPT_STDCALL nlopt_set_xtol_abs(nlopt_opt opt, const double *xtol_abs)
{
     if (opt) {
	  memcpy(opt->xtol_abs, xtol_abs, opt->n * sizeof(double));
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result
NLOPT_STDCALL nlopt_set_xtol_abs1(nlopt_opt opt, double xtol_abs)
{
     if (opt) {
	  unsigned i;
	  for (i = 0; i < opt->n; ++i)
	       opt->xtol_abs[i] = xtol_abs;
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

nlopt_result
NLOPT_STDCALL nlopt_get_xtol_abs(const nlopt_opt opt, double *xtol_abs)
{
     memcpy(xtol_abs, opt->xtol_abs, opt->n * sizeof(double));
     return NLOPT_SUCCESS;
}

GETSET(maxeval, int, maxeval)

/*************************************************************************/

nlopt_result
NLOPT_STDCALL nlopt_set_force_stop(nlopt_opt opt, int force_stop)
{
     if (opt) {
	  opt->force_stop = force_stop;
	  if (opt->force_stop_child)
	       return nlopt_set_force_stop(opt->force_stop_child, force_stop);
	  return NLOPT_SUCCESS;
     }
     return NLOPT_INVALID_ARGS;
}

GET(force_stop, int, force_stop)
nlopt_result NLOPT_STDCALL nlopt_force_stop(nlopt_opt opt) { 
     return nlopt_set_force_stop(opt, 1); 
}

/*************************************************************************/

GET(algorithm, nlopt_algorithm, algorithm)
GET(dimension, unsigned, n)

/*************************************************************************/

GETSET(vector_storage, unsigned, vector_storage)

/*************************************************************************/

nlopt_result NLOPT_STDCALL nlopt_set_initial_step1(nlopt_opt opt, double dx)
{
     unsigned i;
     if (!opt || dx == 0) return NLOPT_INVALID_ARGS;
     if (!opt->dx && opt->n > 0) {
	  opt->dx = (double *) malloc(sizeof(double) * (opt->n));
	  if (!opt->dx) return NLOPT_OUT_OF_MEMORY;
     }
     for (i = 0; i < opt->n; ++i) opt->dx[i] = dx;
     return NLOPT_SUCCESS;
}

nlopt_result
NLOPT_STDCALL nlopt_set_initial_step(nlopt_opt opt, const double *dx)
{
     unsigned i;
     if (!opt) return NLOPT_INVALID_ARGS;
     if (!dx) {
	  free(opt->dx); opt->dx = NULL;
	  return NLOPT_SUCCESS;
     }
     for (i = 0; i < opt->n; ++i) if (dx[i] == 0) return NLOPT_INVALID_ARGS;
     if (!opt->dx && nlopt_set_initial_step1(opt, 1) == NLOPT_OUT_OF_MEMORY)
          return NLOPT_OUT_OF_MEMORY;
     memcpy(opt->dx, dx, sizeof(double) * (opt->n));
     return NLOPT_SUCCESS;
}

nlopt_result
NLOPT_STDCALL nlopt_set_default_initial_step(nlopt_opt opt, const double *x)
{
     const double *lb, *ub;
     unsigned i;

     if (!opt || !x) return NLOPT_INVALID_ARGS;
     lb = opt->lb; ub = opt->ub;

     if (!opt->dx && nlopt_set_initial_step1(opt, 1) == NLOPT_OUT_OF_MEMORY)
	  return NLOPT_OUT_OF_MEMORY;

     /* crude heuristics for initial step size of nonderivative algorithms */
     for (i = 0; i < opt->n; ++i) {
	  double step = HUGE_VAL;

	  if (!nlopt_isinf(ub[i]) && !nlopt_isinf(lb[i])
	      && (ub[i] - lb[i]) * 0.25 < step && ub[i] > lb[i])
	       step = (ub[i] - lb[i]) * 0.25;
	  if (!nlopt_isinf(ub[i]) 
	      && ub[i] - x[i] < step && ub[i] > x[i])
	       step = (ub[i] - x[i]) * 0.75;
	  if (!nlopt_isinf(lb[i]) 
	      && x[i] - lb[i] < step && x[i] > lb[i])
	       step = (x[i] - lb[i]) * 0.75;

	  if (nlopt_isinf(step)) {
	       if (!nlopt_isinf(ub[i]) 
		   && fabs(ub[i] - x[i]) < fabs(step))
		    step = (ub[i] - x[i]) * 1.1;
	       if (!nlopt_isinf(lb[i]) 
		   && fabs(x[i] - lb[i]) < fabs(step))
		    step = (x[i] - lb[i]) * 1.1;
	  }
	  if (nlopt_isinf(step) || step == 0) {
	       step = x[i];
	  }
	  if (nlopt_isinf(step) || step == 0)
	       step = 1;
	  
	  opt->dx[i] = step;
     }
     return NLOPT_SUCCESS;
}

/*************************************************************************/

void NLOPT_STDCALL nlopt_set_munge(nlopt_opt opt,
				   nlopt_munge munge_on_destroy,
				   nlopt_munge munge_on_copy) {
     if (opt) {
	  opt->munge_on_destroy = munge_on_destroy;
	  opt->munge_on_copy = munge_on_copy;
     }
}

void NLOPT_STDCALL nlopt_munge_data(nlopt_opt opt,
                                    nlopt_munge2 munge, void *data) {
     if (opt && munge) {
          unsigned i;
          opt->f_data = munge(opt->f_data, data);
          for (i = 0; i < opt->m; ++i)
               opt->fc[i].f_data = munge(opt->fc[i].f_data, data);
          for (i = 0; i < opt->p; ++i)
               opt->h[i].f_data = munge(opt->h[i].f_data, data);
     }
}



#include <math.h>
#include <float.h>

/*************************************************************************/

int nlopt_isinf(double x) {
    return fabs(x) >= HUGE_VAL * 0.99
#ifdef HAVE_ISINF
    || isinf(x)
#endif
    ;
}
/*************************************************************************/

/*************************************************************************/

static const char nlopt_algorithm_names[NLOPT_NUM_ALGORITHMS][256] = {
    "COBYLA (Constrained Optimization BY Linear Approximations) (local, no-derivative)"
};

const char * NLOPT_STDCALL nlopt_algorithm_name(nlopt_algorithm a)
{
    return nlopt_algorithm_names[a];
}




/* Return a new array of length n (> 0) that gives a rescaling factor
 for each dimension, or NULL if out of memory, with dx being the
 array of nonzero initial steps in each dimension.  */
double *nlopt_compute_rescaling(unsigned n, const double *dx)
{
    double *s = (double *) malloc(sizeof(double) * n);
    unsigned i;
    
    if (!s) return NULL;
    for (i = 0; i < n; ++i) s[i] = 1.0; /* default: no rescaling */
    if (n == 1) return s;
    
    for (i = 1; i < n && dx[i] == dx[i-1]; ++i) ;
    if (i < n) { /* unequal initial steps, rescale to make equal to dx[0] */
        for (i = 1; i < n; ++i)
            s[i] = dx[i] / dx[0];
    }
    return s;
}

void nlopt_rescale(unsigned n, const double *s, const double *x, double *xs)
{
    unsigned i;
    if (!s) { for (i = 0; i < n;++i) xs[i] = x[i]; }
    else { for (i = 0; i < n;++i) xs[i] = x[i] / s[i]; }
}

void nlopt_unscale(unsigned n, const double *s, const double *x, double *xs)
{
    unsigned i;
    if (!s) { for (i = 0; i < n;++i) xs[i] = x[i]; }
    else { for (i = 0; i < n;++i) xs[i] = x[i] * s[i]; }
}

/* return a new array of length n equal to the original array
 x divided by the scale factors s, or NULL on a memory error */
double *nlopt_new_rescaled(unsigned n, const double *s, const double *x)
{
    double *xs = (double *) malloc(sizeof(double) * n);
    if (!xs) return NULL;
    nlopt_rescale(n, s, x, xs);
    return xs;
}

/* since rescaling can flip the signs of the x components and the bounds,
 we may have to re-order the bounds in order to ensure that they
 remain in the correct order */
void nlopt_reorder_bounds(unsigned n, double *lb, double *ub)
{
    unsigned i;
    for (i = 0; i < n; ++i)
        if (lb[i] > ub[i]) {
            double t = lb[i];
            lb[i] = ub[i];
            ub[i] = t;
        }
}


/* utility routines to implement the various stopping criteria */

static int relstop(double vold, double vnew, double reltol, double abstol)
{
    if (nlopt_isinf(vold)) return 0;
    return(fabs(vnew - vold) < abstol
           || fabs(vnew - vold) < reltol * (fabs(vnew) + fabs(vold)) * 0.5
           || (reltol > 0 && vnew == vold)); /* catch vnew == vold == 0 */
}

int nlopt_stop_ftol(const nlopt_stopping *s, double f, double oldf)
{
    return (relstop(oldf, f, s->ftol_rel, s->ftol_abs));
}

int nlopt_stop_f(const nlopt_stopping *s, double f, double oldf)
{
    return (f <= s->minf_max || nlopt_stop_ftol(s, f, oldf));
}

int nlopt_stop_x(const nlopt_stopping *s, const double *x, const double *oldx)
{
    unsigned i;
    for (i = 0; i < s->n; ++i)
        if (!relstop(oldx[i], x[i], s->xtol_rel, s->xtol_abs[i]))
            return 0;
    return 1;
}

int nlopt_stop_dx(const nlopt_stopping *s, const double *x, const double *dx)
{
    unsigned i;
    for (i = 0; i < s->n; ++i)
        if (!relstop(x[i] - dx[i], x[i], s->xtol_rel, s->xtol_abs[i]))
            return 0;
    return 1;
}

static double sc(double x, double smin, double smax)
{
    return smin + x * (smax - smin);
}

/* some of the algorithms rescale x to a unit hypercube, so we need to
 scale back before we can compare to the tolerances */
int nlopt_stop_xs(const nlopt_stopping *s,
                  const double *xs, const double *oldxs,
                  const double *scale_min, const double *scale_max)
{
    unsigned i;
    for (i = 0; i < s->n; ++i)
        if (relstop(sc(oldxs[i], scale_min[i], scale_max[i]),
                    sc(xs[i], scale_min[i], scale_max[i]),
                    s->xtol_rel, s->xtol_abs[i]))
            return 1;
    return 0;
}

int nlopt_stop_evals(const nlopt_stopping *s)
{
    return (s->maxeval > 0 && s->nevals >= s->maxeval);
}




int nlopt_stop_forced(const nlopt_stopping *stop)
{
    return stop->force_stop && *(stop->force_stop);
}

unsigned nlopt_count_constraints(unsigned p, const nlopt_constraint *c)
{
    unsigned i, count = 0;
    for (i = 0; i < p; ++i)
        count += c[i].m;
    return count;
}

unsigned nlopt_max_constraint_dim(unsigned p, const nlopt_constraint *c)
{
    unsigned i, max_dim = 0;
    for (i = 0; i < p; ++i)
        if (c[i].m > max_dim)
            max_dim = c[i].m;
    return max_dim;
}

void nlopt_eval_constraint(double *result, double *grad,
                           const nlopt_constraint *c,
                           unsigned n, const double *x)
{
    if (c->f)
        result[0] = c->f(n, x, grad, c->f_data);
    else
        c->mf(c->m, result, n, x, grad, c->f_data);
}



/*********************************************************************/


/*********************************************************************/
/* wrapper functions, only for derivative-free methods, that
 eliminate dimensions with lb == ub.   (The gradient-based methods
 should handle this case directly, since they operate on much
 larger vectors where I am loathe to make copies unnecessarily.) */

typedef struct {
    nlopt_func f;
    nlopt_mfunc mf;
    void *f_data;
    unsigned n; /* true dimension */
    double *x; /* scratch vector of length n */
    double *grad; /* optional scratch vector of length n */
    const double *lb, *ub; /* bounds, of length n */
} elimdim_data;





/* compute the eliminated dimension: number of dims with lb[i] != ub[i] */
static unsigned elimdim_dimension(unsigned n, const double *lb, const double *ub)
{
    unsigned n0 = 0, i;
    for (i = 0; i < n; ++i) n0 += lb[i] != ub[i] ? 1U : 0;
    return n0;
}



/* inverse of elimdim_shrink */
static void elimdim_expand(unsigned n, double *v,
                           const double *lb, const double *ub)
{
    unsigned i, j;
    if (v && n > 0) {
        j = elimdim_dimension(n, lb, ub) - 1;
        for (i = n - 1; i > 0; --i) {
            if (lb[i] != ub[i])
                v[i] = v[j--];
            else
                v[i] = lb[i];
        }
        if (lb[0] == ub[0])
            v[0] = lb[0];
    }
}


/* like nlopt_destroy, but also frees elimdim_data */
static void elimdim_destroy(nlopt_opt opt)
{
    unsigned i;
    if (!opt) return;
    
    free(((elimdim_data*) opt->f_data)->x);
    free(((elimdim_data*) opt->f_data)->grad);
    free(opt->f_data); opt->f_data = NULL;
    
    for (i = 0; i < opt->m; ++i) {
        free(opt->fc[i].f_data);
        opt->fc[i].f_data = NULL;
    }
    for (i = 0; i < opt->p; ++i) {
        free(opt->h[i].f_data);
        opt->h[i].f_data = NULL;
    }
    
    nlopt_destroy(opt);
}



/*********************************************************************/

#define POP(defaultpop) (opt->stochastic_population > 0 ?		\
opt->stochastic_population :			\
(nlopt_stochastic_population > 0 ?		\
nlopt_stochastic_population : (defaultpop)))

/* unlike nlopt_optimize() below, only handles minimization case */
static nlopt_result nlopt_optimize_(nlopt_opt opt, double *x, double *minf)
{
    const double *lb, *ub;
    nlopt_algorithm algorithm;
    nlopt_func f; void *f_data;
    unsigned n, i;
    nlopt_stopping stop;
    
    if (!opt || !x || !minf || !opt->f
        || opt->maximize) return NLOPT_INVALID_ARGS;
    
    /* reset stopping flag */
    nlopt_set_force_stop(opt, 0);
    opt->force_stop_child = NULL;
    
    /* copy a few params to local vars for convenience */
    n = opt->n;
    lb = opt->lb; ub = opt->ub;
    algorithm = opt->algorithm;
    f = opt->f; f_data = opt->f_data;
    
    if (n == 0) { /* trivial case: no degrees of freedom */
        *minf = opt->f(n, x, NULL, opt->f_data);
        return NLOPT_SUCCESS;
    }
    
    *minf = HUGE_VAL;
    
    
    /* check bound constraints */
    for (i = 0; i < n; ++i)
        if (lb[i] > ub[i] || x[i] < lb[i] || x[i] > ub[i])
            return NLOPT_INVALID_ARGS;
    
    stop.n = n;
    stop.minf_max = opt->stopval;
    stop.ftol_rel = opt->ftol_rel;
    stop.ftol_abs = opt->ftol_abs;
    stop.xtol_rel = opt->xtol_rel;
    stop.xtol_abs = opt->xtol_abs;
    stop.nevals = 0;
    stop.maxeval = opt->maxeval;
    stop.force_stop = &(opt->force_stop);
    
    switch (algorithm) {
            
        case NLOPT_LN_COBYLA: {
            nlopt_result ret;
            int freedx = 0;
            if (!opt->dx) {
                if (nlopt_set_default_initial_step(opt, x) != NLOPT_SUCCESS)
                    return NLOPT_OUT_OF_MEMORY;
            }
            return cobyla_minimize(n, f, f_data,
                                   opt->m, opt->fc,
                                   opt->p, opt->h,
                                   lb, ub, x, minf, &stop,
                                   opt->dx);
            if (freedx) { free(opt->dx); opt->dx = NULL; }
            return ret;
        }
            
            
        default:
            return NLOPT_INVALID_ARGS;
    }
    
    return NLOPT_SUCCESS; /* never reached */
}

/*********************************************************************/

typedef struct {
    nlopt_func f;
    nlopt_precond pre;
    void *f_data;
} f_max_data;

/* wrapper for maximizing: just flip the sign of f and grad */
static double f_max(unsigned n, const double *x, double *grad, void *data)
{
    f_max_data *d = (f_max_data *) data;
    double val = d->f(n, x, grad, d->f_data);
    if (grad) {
        unsigned i;
        for (i = 0; i < n; ++i)
            grad[i] = -grad[i];
    }
    return -val;
}

static void pre_max(unsigned n, const double *x, const double *v,
                    double *vpre, void *data)
{
    f_max_data *d = (f_max_data *) data;
    unsigned i;
    d->pre(n, x, v, vpre, d->f_data);
    for (i = 0; i < n; ++i) vpre[i] = -vpre[i];
}

nlopt_result
NLOPT_STDCALL nlopt_optimize(nlopt_opt opt, double *x, double *opt_f)
{
    nlopt_func f; void *f_data; nlopt_precond pre;
    f_max_data fmd;
    int maximize;
    nlopt_result ret;
    
    if (!opt || !opt_f || !opt->f) return NLOPT_INVALID_ARGS;
    f = opt->f; f_data = opt->f_data; pre = opt->pre;
    
    /* for maximizing, just minimize the f_max wrapper, which
     flips the sign of everything */
    if ((maximize = opt->maximize)) {
        fmd.f = f; fmd.f_data = f_data; fmd.pre = pre;
        opt->f = f_max; opt->f_data = &fmd;
        if (opt->pre) opt->pre = pre_max;
        opt->stopval = -opt->stopval;
        opt->maximize = 0;
    }
    
    { /* possibly eliminate lb == ub dimensions for some algorithms */
        nlopt_opt elim_opt = opt;
        
        ret = nlopt_optimize_(elim_opt, x, opt_f);
        
        if (elim_opt != opt) {
            elimdim_destroy(elim_opt);
            elimdim_expand(opt->n, x, opt->lb, opt->ub);
        }
    }
    
    if (maximize) { /* restore original signs */
        opt->maximize = maximize;
        opt->stopval = -opt->stopval;
        opt->f = f; opt->f_data = f_data; opt->pre = pre;
     	  *opt_f = -*opt_f;
    }
    
    return ret;
}

/*********************************************************************/

nlopt_result nlopt_optimize_limited(nlopt_opt opt, double *x, double *minf,
                                    int maxeval, double maxtime)
{
    int save_maxeval;
    nlopt_result ret;
    
    if (!opt) return NLOPT_INVALID_ARGS;
    
    save_maxeval = nlopt_get_maxeval(opt);
    
    /* override opt limits if maxeval and/or maxtime are more stringent */
    if (save_maxeval <= 0 || (maxeval > 0 && maxeval < save_maxeval))
        nlopt_set_maxeval(opt, maxeval);
    
    ret = nlopt_optimize(opt, x, minf);
    
    nlopt_set_maxeval(opt, save_maxeval);
    
    return ret;
}

/*********************************************************************/

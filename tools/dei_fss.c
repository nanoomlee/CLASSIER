#include "dei_fss.h"

/**
 * Initialize the integrator_fss
 *
 */
int initialize_generic_integrator_fss(
				  int n_dim,
				  struct generic_integrator_fss_workspace * pgi){

  /** - Allocate workspace dynamically */

  pgi->n = n_dim;

  class_alloc(pgi->y,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->dydx,
	      sizeof(double)*n_dim,
	      pgi->error_message);

  return _SUCCESS_;
}


/**
 * Free the integrator_fss's memory space
 *
 * Called by background_solve(); thermodynamics_solve_with_recfast(); perturbations_solve().
 */
int cleanup_generic_integrator_fss(struct generic_integrator_fss_workspace * pgi){

  free(pgi->y);
  free(pgi->dydx);

  return _SUCCESS_;
}

int generic_integrator_fss(int (*derivs)(double x, double y[], double yprime[], void * parameters_and_workspace, ErrorMsg error_message),
		       double x1,
		       double x2,
		       double ystart[],
		       void * parameters_and_workspace_for_derivs,
		       double eps,
		       double hmin,
		       struct generic_integrator_fss_workspace * pgi)

{
  int i;
  double h = x2-x1;
  class_call((*derivs)(x1,pgi->y,pgi->dydx,parameters_and_workspace_for_derivs, pgi->error_message),
       pgi->error_message,
       pgi->error_message);
  for (i=0;i<pgi->n;i++) ystart[i] += h * pgi->dydx[i];
}


#ifndef __DEI__
#define __DEI__

#include "common.h"

struct generic_integrator_fss_workspace
{
  int n;
  double * y;
  double * dydx;

  /**
    * zone for writing error messages
    */
  ErrorMsg error_message;

};

/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int initialize_generic_integrator_fss(
				    int n_dim,
				    struct generic_integrator_fss_workspace * pgi
				    );

  int cleanup_generic_integrator_fss(struct generic_integrator_fss_workspace * pgi);

  int generic_integrator_fss(int (*derivs)(double x,
				       double y[],
				       double yprime[],
				       void * parameters_and_workspace,
				       ErrorMsg error_message),
			 double x1,
			 double x2,
			 double ystart[],
			 void * parameters_and_workspace_for_derivs,
			 double eps,
			 double hmin,
			 struct generic_integrator_fss_workspace * pgi);

#ifdef __cplusplus
}
#endif

#endif

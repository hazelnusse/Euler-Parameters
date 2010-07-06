/*
 * =====================================================================================
 *
 *       Filename:  rigidbodyeoms.h
 *
 *    Description:  function prototypes defining the system of first order ordinary
 *    differential equations which define the motion of a rigid body with
 *    arbitrarily applied torques and forces.  Makes use of Euler parameters to
 *    define the orientation.
 *
 *         Author:  Dale Lukas Peterson
 *        Company:  University of California, Davis
 *
 * =====================================================================================
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

typedef struct {
  double g, ma, Ixx, Iyy, Izz, Ixy, Iyz, Ixz;
  double Tax, Tay, Taz;
  double ke, pe, te;
  double x[7], f[7], z[47];
  // 4x4 transformation matrix
  double m[16];
  // A and B matrices
  double A[49], B[21];
  // Frame rate to animate at, also controls output points of numerical
  // integration
  double t, tf, h, fps;
  int k, status;
  
  // For storing filename roots passed through command line
  char *pngs;

  // Boiler plate code to use GSL ODE integrator
  const gsl_odeiv_step_type * T;
  gsl_odeiv_step * s;
  gsl_odeiv_control * c;
  gsl_odeiv_evolve * e;
  gsl_odeiv_system sys;
} RigidBody;

int eoms(const double t, const double *x, double f[], void *params);
void evalOutputs(RigidBody * body);
void initRigidBody(RigidBody * body);
void freeRigidBody(RigidBody * body);
void processOptions(int argc, char ** argv, RigidBody * body);

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

typedef struct {
  double ma, Ixx, Iyy, Izz, Ixz;
  double Fax, Fay, Faz, Fnx, Fny, Fnz, Tax, Tay, Taz, Tnx, Tny, Tnz;
  double state[13], f[13];
  double m[16];
} RigidBody;

int eoms_f(const double t, const double *x, double f[], void *params);
void evalTransformation(RigidBody * body);
void initRigidBody(RigidBody * body);


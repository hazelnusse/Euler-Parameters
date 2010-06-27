/*
 * =====================================================================================
 *
 *       Filename:  rigidbodyeoms.c
 *
 *    Description:  Functions defining the system of first order ordinary
 *    differential equations which define the motion of a rigid body with
 *    arbitrarily applied torques and forces.  Makes use of Euler parameters to
 *    define the orientation.
 *
 *         Author:  Dale Lukas Peterson
 *        Company:  University of California, Davis
 *
 * =====================================================================================
 */

#include "rigidbodyeoms.h"
#include <math.h>
#include <gsl/gsl_errno.h>

int eoms_f(const double t, const double *x, double f[], void *params)
{
  // state ordering:  [e0, e1, e2, e3,  x,  y,  z, wx, wy, wz, vx, vy, vz]
  //                    0   1   2   3   4   5   6   7   8   9  10  11  12
  RigidBody * body = (RigidBody *) params;

  f[0] = 0.5*(x[1]*x[9] + x[3]*x[7] - x[2]*x[8]);
  f[1] = 0.5*(x[2]*x[7] + x[3]*x[8] - x[0]*x[9]);
  f[2] = 0.5*(x[0]*x[8] + x[3]*x[9] - x[1]*x[7]);
  f[3] = -0.5*(x[0]*x[7] + x[1]*x[8] + x[2]*x[9]);
  f[4] = (body->Fax+2*body->Fny*(x[0]*x[1]+x[2]*x[3])+2*body->Fnz*(x[0]*x[2]-x[1]*x[3])+body->Fnx*(1-2*pow(x[1],2)-2*pow(x[2],2)))/body->ma + x[9]*x[11] - x[8]*x[12];
  f[5] = (body->Fay+2*body->Fnz*(x[0]*x[3]+x[1]*x[2])+2*body->Fnx*(x[0]*x[1]-x[2]*x[3])+body->Fny*(1-2*pow(x[0],2)-2*pow(x[2],2)))/body->ma + x[7]*x[12] - x[9]*x[10];
  f[6] = x[8]*x[10] - (2*body->Fny*(x[0]*x[3]-x[1]*x[2])-body->Faz-2*body->Fnx*(x[0]*x[2]+x[1]*x[3])-body->Fnz*(1-2*pow(x[0],2)-2*pow(x[1],2)))/body->ma - x[7]*x[11];
  f[7] = 2*(x[0]*x[2]+x[1]*x[3])*x[12] + 2*(x[0]*x[1]-x[2]*x[3])*x[11] + (1-2*pow(x[1],2)-2*pow(x[2],2))*x[10];
  f[8] = 2*(x[0]*x[1]+x[2]*x[3])*x[10] + (1-2*pow(x[0],2)-2*pow(x[2],2))*x[11] - 2*(x[0]*x[3]-x[1]*x[2])*x[12];
  f[9] = 2*(x[0]*x[3]+x[1]*x[2])*x[11] + 2*(x[0]*x[2]-x[1]*x[3])*x[10] + (1-2*pow(x[0],2)-2*pow(x[1],2))*x[12];
  f[10] = (body->Izz*(body->Tax+2*body->Tny*(x[0]*x[1]+x[2]*x[3])+2*body->Tnz*(x[0]*x[2]-x[1]*x[3])+body->Tnx*(1-2*pow(x[1],2)-2*pow(x[2],2))+x[8]*(body->Iyy*x[9]-body->Ixz*x[7]-body->Izz*x[9]))+body->Ixz*(2*body->Tny*(x[0]*x[3]-x[1]*x[2])-body->Taz-2*body->Tnx*(
  x[0]*x[2]+x[1]*x[3])-body->Tnz*(1-2*pow(x[0],2)-2*pow(x[1],2))-x[8]*(body->Ixx*x[7]+body->Ixz*x[9]-body->Iyy*x[7])))/(body->Ixx*body->Izz-pow(body->Ixz,2));
  f[11] = (body->Tay+2*body->Tnz*(x[0]*x[3]+x[1]*x[2])+2*body->Tnx*(x[0]*x[1]-x[2]*x[3])+body->Tny*(1-2*pow(x[0],2)-2*pow(x[2],2))+x[7]*(body->Ixz*x[7]+body->Izz*x[9])-x[9]*(body->Ixx*x[7]+body->Ixz*x[9]))/body->Iyy;
  f[12] = -(body->Ixz*(body->Tax+2*body->Tny*(x[0]*x[1]+x[2]*x[3])+2*body->Tnz*(x[0]*x[2]-x[1]*x[3])+body->Tnx*(1-2*pow(x[1],2)-2*pow(x[2],2))+x[8]*(body->Iyy*x[9]-body->Ixz*x[7]-body->Izz*x[9]))+body->Ixx*(2*body->Tny*(x[0]*x[3]-x[1]*x[2])-body->Taz-2*body->Tnx*(x[0]*x[2]+x[1]*x[3])-body->Tnz*(1-2*pow(x[0],2)-2*pow(x[1],2))-x[8]*(body->Ixx*x[7]+body->Ixz*x[9]-body->Iyy*x[7])))/(body->Ixx*body->Izz-pow(body->Ixz,2));

  return GSL_SUCCESS;
} // eoms_f()

void evalTransformation(RigidBody * body)
{
  body->m[0] = 1 - 2*pow(body->state[1],2) - 2*pow(body->state[2],2);
  body->m[1] = 2*body->state[0]*body->state[1] + 2*body->state[2]*body->state[3];
  body->m[2] = 2*body->state[0]*body->state[2] - 2*body->state[1]*body->state[3];
  body->m[4] = 2*body->state[0]*body->state[1] - 2*body->state[2]*body->state[3];
  body->m[5] = 1 - 2*pow(body->state[0],2) - 2*pow(body->state[2],2);
  body->m[6] = 2*body->state[0]*body->state[3] + 2*body->state[1]*body->state[2];
  body->m[8] = 2*body->state[0]*body->state[2] + 2*body->state[1]*body->state[3];
  body->m[9] = 2*body->state[1]*body->state[2] - 2*body->state[0]*body->state[3];
  body->m[10] = 1 - 2*pow(body->state[0],2) - 2*pow(body->state[1],2);
  body->m[12] = body->state[4];
  body->m[13] = body->state[5];
  body->m[14] = body->state[6];
} // evalTransformation()

void initRigidBody(RigidBody * body)
{
  body->m[3] = body->m[7] = body->m[11] = 0.0;
  body->m[15] = 1.0;
} // initRigidBody()

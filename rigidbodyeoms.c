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

#include <math.h>
#include <getopt.h>
#include <string.h>

#include "rigidbodyeoms.h"

int eoms(const double t, const double *VAR, double VARp[], void *params)
{
  // state ordering:  [e0, e1, e2, e3,  x,  y,  z, u0, u1, u2, u3, u4, u5]
  //                    0   1   2   3   4   5   6   7   8   9  10  11  12
  RigidBody * body = (RigidBody *) params;
  double  e0, e1, e2, e3, u0, u1, u2, e0p, e1p, e2p, e3p, u0p, u1p, u2p;

  double Ixx = body->Ixx, Iyy = body->Iyy, Izz = body->Izz, Ixy = body->Ixy, Iyz = body->Iyz, Ixz = body->Ixz;
  double Tax = body->Tax, Tay = body->Tay, Taz = body->Taz;

  // Determine magnitude of quaternion
  double mag = 1.0; sqrt(VAR[0]*VAR[0] + VAR[1]*VAR[1] + VAR[2]*VAR[2] + VAR[3]*VAR[3]);
  double * z = body->z;

  /* Evaluate constants */
  z[21] = Ixy*Iyz - Ixz*Iyy;
  z[20] = Ixx*Iyz - Ixy*Ixz;
  z[19] = Ixx*Iyy - pow(Ixy,2);
  z[22] = Iyz*z[20] - Ixz*z[21] - Izz*z[19];
  z[23] = Iyy*Izz - pow(Iyz,2);
  z[24] = Ixy*Izz - Ixz*Iyz;
  z[25] = Ixx*Izz - pow(Ixz,2);

/* Update variables after integration step */
  e0 = VAR[0]/mag;
  e1 = VAR[1]/mag;
  e2 = VAR[2]/mag;
  e3 = VAR[3]/mag;
  u0 = VAR[4];
  u1 = VAR[5];
  u2 = VAR[6];

  e0p = 0.5*e1*u2 + 0.5*e3*u0 - 0.5*e2*u1;
  e1p = 0.5*e2*u0 + 0.5*e3*u1 - 0.5*e0*u2;
  e2p = 0.5*e0*u1 + 0.5*e3*u2 - 0.5*e1*u0;
  e3p = -0.5*e0*u0 - 0.5*e1*u1 - 0.5*e2*u2;
  z[11] = Ixy*u0 + Iyy*u1 + Iyz*u2;
  z[10] = Ixx*u0 + Ixy*u1 + Ixz*u2;
  z[13] = u0*z[11] - u1*z[10];
  z[18] = z[13] - Taz;
  z[12] = Ixz*u0 + Iyz*u1 + Izz*u2;
  z[15] = u1*z[12] - u2*z[11];
  z[16] = z[15] - Tax;
  z[14] = u2*z[10] - u0*z[12];
  z[17] = z[14] - Tay;
  z[26] = (z[21]*z[18]+z[23]*z[16]-z[24]*z[17])/z[22];
  u0p = z[26];
  z[27] = (z[20]*z[18]+z[24]*z[16]-z[25]*z[17])/z[22];
  u1p = -z[27];
  z[28] = (z[20]*z[17]-z[19]*z[18]-z[21]*z[16])/z[22];
  u2p = -z[28];

/* Update derivative array prior to integration step */
  VARp[0] = e0p;
  VARp[1] = e1p;
  VARp[2] = e2p;
  VARp[3] = e3p;
  VARp[4] = u0p;
  VARp[5] = u1p;
  VARp[6] = u2p;
  // Return sucess
  return GSL_SUCCESS;
} // eoms()

void evalOutputs(RigidBody * body)
{
  double  e0 = body->x[0], e1 = body->x[1], e2 = body->x[2], e3 = body->x[3], u0 = body->x[4], u1 = body->x[5], u2 = body->x[6];
  double * m = body->m, * A = body->A, * B = body->B;

  double Ixx = body->Ixx, Iyy = body->Iyy, Izz = body->Izz, Ixy = body->Ixy, Iyz = body->Iyz, Ixz = body->Ixz;
  double * z = body->z;

/* Evaluate output quantities */
  z[1] = 1 - 2*pow(e1,2) - 2*pow(e2,2);
  z[2] = 2*e0*e1 - 2*e2*e3;
  z[3] = 2*e0*e2 + 2*e1*e3;
  z[4] = 2*e0*e1 + 2*e2*e3;
  z[5] = 1 - 2*pow(e0,2) - 2*pow(e2,2);
  z[6] = 2*e1*e2 - 2*e0*e3;
  z[7] = 2*e0*e2 - 2*e1*e3;
  z[8] = 2*e0*e3 + 2*e1*e2;
  z[9] = 1 - 2*pow(e0,2) - 2*pow(e1,2);
  z[29] = Ixy*u0 + z[11] - Ixx*u1;
  z[30] = Ixz*u1 - Ixy*u2;
  z[31] = Ixx*u2 - Ixz*u0 - z[12];
  z[32] = (z[21]*z[29]+z[23]*z[30]-z[24]*z[31])/z[22];
  z[33] = Iyy*u0 - Ixy*u1 - z[10];
  z[34] = Iyz*u1 + z[12] - Iyy*u2;
  z[35] = Ixy*u2 - Iyz*u0;
  z[36] = (z[21]*z[33]+z[23]*z[34]-z[24]*z[35])/z[22];
  z[37] = Iyz*u0 - Ixz*u1;
  z[38] = Izz*u1 - Iyz*u2 - z[11];
  z[39] = Ixz*u2 + z[10] - Izz*u0;
  z[40] = (z[21]*z[37]+z[23]*z[38]-z[24]*z[39])/z[22];
  z[41] = (z[20]*z[29]+z[24]*z[30]-z[25]*z[31])/z[22];
  z[42] = (z[20]*z[33]+z[24]*z[34]-z[25]*z[35])/z[22];
  z[43] = (z[20]*z[37]+z[24]*z[38]-z[25]*z[39])/z[22];
  z[44] = (z[20]*z[31]-z[19]*z[29]-z[21]*z[30])/z[22];
  z[45] = (z[20]*z[35]-z[19]*z[33]-z[21]*z[34])/z[22];
  z[46] = (z[20]*z[39]-z[19]*z[37]-z[21]*z[38])/z[22];

  m[0] = z[1];
  m[1] = z[4];
  m[2] = z[7];
  m[3] = 0;
  m[4] = z[2];
  m[5] = z[5];
  m[6] = z[8];
  m[7] = 0;
  m[8] = z[3];
  m[9] = z[6];
  m[10] = z[9];
  m[11] = 0;
  m[12] = 0;
  m[13] = 0;
  m[14] = 0;
  m[15] = 1;
  A[0] = 0;
  A[1] = 0.5*u2;
  A[2] = -0.5*u1;
  A[3] = 0.5*u0;
  A[4] = 0.5*e3;
  A[5] = -0.5*e2;
  A[6] = 0.5*e1;
  A[7] = -0.5*u2;
  A[8] = 0;
  A[9] = 0.5*u0;
  A[10] = 0.5*u1;
  A[11] = 0.5*e2;
  A[12] = 0.5*e3;
  A[13] = -0.5*e0;
  A[14] = 0.5*u1;
  A[15] = -0.5*u0;
  A[16] = 0;
  A[17] = 0.5*u2;
  A[18] = -0.5*e1;
  A[19] = 0.5*e0;
  A[20] = 0.5*e3;
  A[21] = -0.5*u0;
  A[22] = -0.5*u1;
  A[23] = -0.5*u2;
  A[24] = 0;
  A[25] = -0.5*e0;
  A[26] = -0.5*e1;
  A[27] = -0.5*e2;
  A[28] = 0;
  A[29] = 0;
  A[30] = 0;
  A[31] = 0;
  A[32] = z[32];
  A[33] = z[36];
  A[34] = z[40];
  A[35] = 0;
  A[36] = 0;
  A[37] = 0;
  A[38] = 0;
  A[39] = -z[41];
  A[40] = -z[42];
  A[41] = -z[43];
  A[42] = 0;
  A[43] = 0;
  A[44] = 0;
  A[45] = 0;
  A[46] = -z[44];
  A[47] = -z[45];
  A[48] = -z[46];
  B[0] = 0;
  B[1] = 0;
  B[2] = 0;
  B[3] = 0;
  B[4] = 0;
  B[5] = 0;
  B[6] = 0;
  B[7] = 0;
  B[8] = 0;
  B[9] = 0;
  B[10] = 0;
  B[11] = 0;
  B[12] = -z[23]/z[22];
  B[13] = z[24]/z[22];
  B[14] = -z[21]/z[22];
  B[15] = z[24]/z[22];
  B[16] = -z[25]/z[22];
  B[17] = z[20]/z[22];
  B[18] = -z[21]/z[22];
  B[19] = z[20]/z[22];
  B[20] = -z[19]/z[22];
} // evalOutputs()

void initRigidBody(RigidBody * body)
{
  int i;
  // Initialize all states to zero
  for (i = 0; i < 7; ++i) {
    body->x[i] = 0.0;
  } // for i
  // Last Euler parameter == 1 aligns body frame with inertial frame
  body->x[3] = 1.0;

  // Initialize mass and inertias
  body->ma = 1.0;
  body->Ixx = 1.0;
  body->Iyy = 2.0;
  body->Izz = 1.0;
  body->Ixy = 0.0;
  body->Iyz = 0.0;
  body->Ixz = 0.0;

  // Initialize all forces and moments to zero, and gravitational field
  body->Tax = body->Tay = body->Taz = 0.0;
  body->g = 0.0;

  // Initialize time to zero, final time to 20.0, frame rate
  body->t = 0.0; body->tf = 20.0; body->fps = 60.0;
  body->h = 1e-3;

  // Initialize frame step counter 0
  body->k = 0;
  body->pngs = NULL;  // by default, don't save pngs

  // Boiler plate code to use GSL ODE integrator
  body->T = gsl_odeiv_step_rk8pd;
  body->s = gsl_odeiv_step_alloc(body->T, 7);
  body->c = gsl_odeiv_control_y_new(1e-6, 0.0);
  body->e = gsl_odeiv_evolve_alloc(7);
  body->sys.function = eoms;
  body->sys.jacobian = NULL;
  body->sys.dimension = 7;
  body->sys.params = body;

  // These entries of the 4x4 transformation matrix are constant
  body->m[3] = body->m[7] = body->m[11] = 0.0;
  body->m[15] = 1.0;
} // initRigidBody()

void freeRigidBody(RigidBody * body)
{
  gsl_odeiv_evolve_free(body->e);
  gsl_odeiv_control_free(body->c);
  gsl_odeiv_step_free(body->s);
  free(body);
} // freeRigidBody

void processOptions(int argc, char ** argv, RigidBody * body)
{
  int c, opt_index;
  struct option long_options[] = {
     {"help", no_argument, 0, '?'},
     {"pngs", required_argument, 0, 'p'},
     {"Ixx",  required_argument, 0, 'a'},
     {"Iyy",  required_argument, 0, 'b'},
     {"Izz",  required_argument, 0, 'c'},
     {"Ixy",  required_argument, 0, 'd'},
     {"Iyz",  required_argument, 0, 'e'},
     {"Ixz",  required_argument, 0, 'f'},
     {"wx",  required_argument, 0, 'g'},
     {"wy",  required_argument, 0, 'h'},
     {"wz",  required_argument, 0, 'i'},
     {"tf",  required_argument, 0, 't'},
     {0, 0, 0, 0} };
  while (1) {
    opt_index = 0;
    c = getopt_long(argc, argv, "?a:b:c:d:e:f:g:h:i:t:p:", long_options, &opt_index);

  if (c == -1)
    break;

  switch (c) {
    case '?':
      printf("usage: %s [OPTION]\n\n"
             "Mandatory arguments to long options are mandatory for short options too.\n\n"
             "  -?, --help                         display this help and exit.\n"
             "  -a val, --Ixx=val                      Specify Ixx moment of inertia.\n"
             "  -b val, --Iyy=val                      Specify Iyy moment of inertia.\n"
             "  -c val, --Izz=val                      Specify Izz moment of inertia.\n"
             "  -d val, --Ixy=val                      Specify Ixy moment of inertia.\n"
             "  -e val, --Iyz=val                      Specify Iyz moment of inertia.\n"
             "  -f val, --Ixz=val                      Specify Ixz moment of inertia.\n"
             "  -g val, --wx=val                       Specify initial angular velocity about body-fixed x axis\n"
             "  -h val, --wy=val                       Specify initial angular velocity about body-fixed y axis\n"
             "  -i val, --wz=val                       Specify initial angular velocity about body-fixed z axis\n"
             "  -t val, --tf=val                       Specify total simulation time\n\n"
             "Example of how to specify Ixx=1.0, Iyy=2.0, Izz=3.0, intial angular velocity of w=[0.1, 2.0, 0.1]:\n\n"
             "$ %s -a 1.0 --Iyy=2.0 --Izz=3.0 -g 0.1 --wy=2.0 --wz=0.1\n\n", 
             argv[0], argv[0]);
      exit(0);
 
    case 'a': body->Ixx = atof(optarg); break;
    case 'b': body->Iyy = atof(optarg); break;
    case 'c': body->Izz = atof(optarg); break;
    case 'd': body->Ixy = atof(optarg); break;
    case 'e': body->Iyz = atof(optarg); break;
    case 'f': body->Ixz = atof(optarg); break;
    case 'g': body->x[4] = atof(optarg); break;
    case 'h': body->x[5] = atof(optarg); break;
    case 'i': body->x[6] = atof(optarg); break;
    case 't': body->tf = atof(optarg); break;
    case 'p': 
      body->pngs = (char *) malloc(strlen(optarg) + 9);
      strcpy(body->pngs, optarg);
      break;
    default: abort();
    } // switch(c)
  } // while
} // processOptions()

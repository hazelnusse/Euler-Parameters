/*
 * =====================================================================================
 *
 *       Filename:  simulate.c
 *
 *    Description:  Example of using the GNU Scientific Library to numerically
 *                  integrate equations of motion of a single rigid body
 *
 *         Author:  Dale Lukas Peterson
 *        Company:  University of California Davis
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include <gsl/gsl_errno.h>
#include "rigidbodyeoms.h"
#define WIDTH 1280
#define HEIGHT 720

// Declare a global pointer to a RigidBody structure
RigidBody * body;

void init(void)
{
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glDisable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT |GL_DEPTH_BUFFER_BIT);
  // glMatrixMode(GL_MODELVIEW);
  
  glLoadIdentity();
  glTranslatef(0.0, 0.0, -3.0);

  //Add ambient light
  GLfloat ambientColor[] = {0.4f, 0.4f, 0.4f, 1.0f}; //Color (0.2, 0.2, 0.2)
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);

  //Add positioned light
  GLfloat lightColor0[] = {0.7f, 0.7f, 0.7f, 1.0f}; //Color (0.5, 0.5, 0.5)
  GLfloat lightPos0[] = {0.0f, 0.0f, -2.0f, 1.0f}; //Positioned at (4, 0, 8)
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);

  //Add directed light
  GLfloat lightColor1[] = {0.5f, 0.2f, 0.2f, 1.0f}; //Color (0.5, 0.2, 0.2)
  //Coming from the direction (-1, 0.5, 0.5)
  GLfloat lightPos1[] = {-1.0f, 0.0f, 0.0f, 0.0f};
  glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);

  glMultMatrixd(body->m);
  // x axis
  /*
  float fvViewMatrix[16]; 
  int i;
  glGetFloatv(GL_MODELVIEW_MATRIX, fvViewMatrix);
  printf("[");
  for (i = 0; i < 3; ++i)
    printf("[%f, %f, %f, %f];\n", body->m[i], body->m[i + 4], body->m[i + 8], body->m[i + 12]);
  printf("[%f, %f, %f, %f]]\n", body->m[3], body->m[7], body->m[11], body->m[15]);
  */
  glBegin(GL_LINES);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
  glEnd();

  glPushMatrix();
    glTranslatef(0.9, 0.0, 0.0);
    glRotatef(90.0, 0.0, 1.0, 0.0);
    glutSolidCone(0.05, 0.1, 10, 1);
  glPopMatrix();

  // y axis
  glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
  glEnd();

  glPushMatrix();
    glTranslatef(0.0, 0.9, 0.0);
    glRotatef(-90.0, 1.0, 0.0, 0.0);
    glutSolidCone(0.05, 0.1, 10, 1);
  glPopMatrix();

  // z axis
  glBegin(GL_LINES);
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
  glEnd();

  glPushMatrix();
    glTranslatef(0.0, 0.0, 0.9);
    glutSolidCone(0.05, 0.1, 10, 1);
  glPopMatrix();

  glutSwapBuffers();  // Only needed if in double buffer mode
}

void updateState(int value)
{
  double tj, mag;

  ++(body->k);
  tj = body->t + (1.0 / body->fps);

  while (body->t < tj) {
    body->status = gsl_odeiv_evolve_apply(body->e, body->c, body->s, &(body->sys), &(body->t), tj, &(body->h), body->x);
    mag = sqrt(body->x[0]*body->x[0] + body->x[1]*body->x[1] + body->x[2]*body->x[2] + body->x[3]*body->x[3]);
    body->x[0] /= mag;
    body->x[1] /= mag;
    body->x[2] /= mag;
    body->x[3] /= mag;
  }

  // Evaluate output quantities
  evalOutputs(body);
  // Print the magnitude of the quaternion
  printf("sqrt(e0^2 + e1^2 + e2^2 + e3^2) = %0.16f\n", sqrt(pow(body->x[0], 2.0) + pow(body->x[1], 2.0) + pow(body->x[2], 2.0) + pow(body->x[3], 2.0)));

  glutPostRedisplay();
  
  if (body->k == floor(body->fps * body->tf))
    exit(0);
  // re-register the callback
  glutTimerFunc((unsigned int) (1000.0/body->fps), updateState, 0);
}

void reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective(45.0, (GLfloat) w/(GLfloat) h, 1.0, 20.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void keyboard (unsigned char key, int x, int y)
{
  if (key == 27) {
    free(body);
    exit(0);
  }
}


int main(int argc, char ** argv)
{
  // Dynamical allocate memory for a RigidBody structure
  body = (RigidBody *) malloc(sizeof(RigidBody));

  initRigidBody(body);  // set some default mass, inertia, forces, initial conditions

  body->Ixx = 1.0;
  body->Iyy = 2.0;
  body->Izz = 3.0;

  // Initial body fixed angular rates
  body->x[4] = 0.1;
  body->x[5] = 3.0;
  body->x[6] = 0.1;
  body->tf = 10.0;
  evalOutputs(body);

  // Initialize animation window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(WIDTH, HEIGHT);
  glutCreateWindow("Euler parameter animation");
  init();

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutTimerFunc((unsigned int) (1000.0/(body->fps)), updateState, 0);
  glutMainLoop();
  return 0;
}

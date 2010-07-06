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

#ifdef __APPLE__                                                                                                                                                                                                
#include <OpenGL/OpenGL.h>                                                                                                                                                                                      
#include <GLUT/glut.h>                                                                                                                                                                                          
#else                                                                                                                                                                                                           
#include <GL/glut.h>                                                                                                                                                                                            
#endif    

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>

#include "rigidbodyeoms.h"
#define WIDTH 720
#define HEIGHT 720

// Declare a global pointer to a RigidBody structure
RigidBody * body;
char wx[15], wy[15], wz[15], t[15];

void render_string( char* string, float x, float y )
{
  int i;
  //Set initial string raster position.  Assumes x-y plane is the plane of the monitor with origin in upper-left corner (I think)
  glRasterPos3f(x, y, 0.0);
  
  //Print string; Each character rendered before raster position iterated by amount required for correct font spacing
  for (i = 0; i < 13; ++i)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *(string+i));
}

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
  render_string( t, -1.0, 1.0);
  render_string(wx, -1.0, 0.9);
  render_string(wy, -1.0, 0.8);
  render_string(wz, -1.0, 0.7);

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
  sprintf( t, " t = %6.1f", body->t);
  sprintf(wx, "wx = %6.5f", body->x[4]);
  sprintf(wy, "wy = %6.5f", body->x[5]);
  sprintf(wz, "wz = %6.5f", body->x[6]);

  // Evaluate output quantities
  evalOutputs(body);
  // Print the magnitude of the quaternion
  // printf("sqrt(e0^2 + e1^2 + e2^2 + e3^2) = %0.16f\n", sqrt(pow(body->x[0], 2.0) + pow(body->x[1], 2.0) + pow(body->x[2], 2.0) + pow(body->x[3], 2.0)));

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
    freeRigidBody(body);
    exit(0);
  }
}

int main(int argc, char ** argv)
{
  // Dynamical allocate memory for a RigidBody structure
  body = (RigidBody *) malloc(sizeof(RigidBody));

  initRigidBody(body);  // set some default mass, inertia, forces, initial conditions

  processOptions(argc, argv, body);
  
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

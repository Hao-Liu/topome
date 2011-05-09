#ifndef _tpm_graphic_h
#define _tpm_graphic_h

#include <sys/time.h>
#include <GL/glx.h>
#include <GL/glu.h>
#include <X11/extensions/xf86vmode.h>
#include "force.h"

#include "types.h"
#include "arcball.h"



//GLfloat rotTri, rotQuad;

void resizeGLScene(unsigned int width, unsigned int height);
int initGL(GLWindow *GLWin);
GLvoid killGLWindow(GLWindow *GLWin);
Bool createGLWindow(char* title, int width, int height, int bits, 
										Bool fullscreenflag, GLWindow *GLWin);
void GraphicOutput(System *system, GLWindow *GLWin);
void 
init_graphic(System *system, GLWindow *gl_window);
void
release_graphic (GLWindow *gl_window);
void graphic_output (System *system, GLWindow *GLWin);
#endif //_tpm_graphic_h

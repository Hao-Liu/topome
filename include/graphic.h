#ifndef _tpm_graphic_h
#define _tpm_graphic_h

#include <sys/time.h>
#include <GL/glx.h>
#include <GL/glu.h>
#include <X11/extensions/xf86vmode.h>
#include "molecule.h"

#include "system.h"
#include "arcball.h"

/* stuff about our window grouped together */
typedef struct {
    Display *dpy;
    int screen;
    Window win;
    GLXContext ctx;
    XSetWindowAttributes attr;
    Bool fs;
    Bool doubleBuffered;
    XF86VidModeModeInfo deskMode;
    int x, y;
    
    int isClicked, isRClicked, isDragging;
    Matrix3fT   LastRot, ThisRot;
    Matrix4fT   Transform;
    Point2fT    MousePt;
    ArcBall			arcball;
    
    unsigned int width, height;
    unsigned int depth;    
		struct timeval lastTickCount;				// Tick Counter
} GLWindow;


//GLfloat rotTri, rotQuad;

void resizeGLScene(unsigned int width, unsigned int height);
int initGL(GLWindow *GLWin);
GLvoid killGLWindow(GLWindow *GLWin);
Bool createGLWindow(char* title, int width, int height, int bits, 
										Bool fullscreenflag, GLWindow *GLWin);
void GraphicOutput(System *system, GLWindow *GLWin);
#endif //_tpm_graphic_h

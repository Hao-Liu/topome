#include "molecule.h"
#include <stdlib.h>
#include <GL/glx.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <X11/extensions/xf86vmode.h>
#include <X11/keysym.h>
#include <arcball.h>
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

/* attributes for a single buffered visual in RGBA format with at least
 * 4 bits per color and a 16 bit depth buffer */
static int attrListSgl[] = {GLX_RGBA, GLX_RED_SIZE, 4, 
    GLX_GREEN_SIZE, 4, 
    GLX_BLUE_SIZE, 4, 
    GLX_DEPTH_SIZE, 16,
    None};

/* attributes for a double buffered visual in RGBA format with at least 
 * 4 bits per color and a 16 bit depth buffer */
static int attrListDbl[] = { GLX_RGBA, GLX_DOUBLEBUFFER, 
    GLX_RED_SIZE, 4, 
    GLX_GREEN_SIZE, 4, 
    GLX_BLUE_SIZE, 4, 
    GLX_DEPTH_SIZE, 16,
    None };

GLWindow GLWin;
GLfloat rotTri, rotQuad;

void resizeGLScene(unsigned int width, unsigned int height);
int initGL(void);
GLvoid killGLWindow(void);
Bool createGLWindow(char* title, int width, int height, int bits, 
										Bool fullscreenflag);
int RenderMolecules(SYSTEM *system);


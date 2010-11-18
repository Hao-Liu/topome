#include "graphic.h"
	float AtomColor[][3]=
	{
	{1.0f, 1.0f, 1.0f}, //X
	{1.0f, 1.0f, 1.0f}, //H
	{1.0f, 1.0f, 1.0f}, //He
	{1.0f, 1.0f, 1.0f}, //Li
	{1.0f, 1.0f, 1.0f}, //Be
	{1.0f, 1.0f, 1.0f}, //B
	{0.5f, 0.5f, 0.5f}, //C
	{0.0f, 0.0f, 1.0f}, //N
	{1.0f, 0.0f, 0.0f}, //O
	{1.0f, 1.0f, 1.0f}, //F
	{1.0f, 1.0f, 1.0f}, //Ne
	{1.0f, 1.0f, 1.0f}, //Na
	{1.0f, 1.0f, 1.0f}, //Mg
	{1.0f, 1.0f, 1.0f}, //Al
	{1.0f, 1.0f, 1.0f}, //Si
	{1.0f, 1.0f, 1.0f}, //P
	{1.0f, 1.0f, 1.0f}, //S
	{1.0f, 1.0f, 1.0f}, //Cl
	{1.0f, 1.0f, 1.0f}, //Ar
	{1.0f, 1.0f, 1.0f}, //K
	{0.5f, 0.5f, 1.0f}, //Ca
	};

void RenderCell(SYSTEM *system)
{
	//Draw cell
	glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(0.0f,0.0f,system->dimension);
		glVertex3f(0.0f,system->dimension,0.0f);
		glVertex3f(0.0f,system->dimension,system->dimension);
		glVertex3f(system->dimension,0.0f,0.0f);
		glVertex3f(system->dimension,0.0f,system->dimension);
		glVertex3f(system->dimension,system->dimension,0.0f);
		glVertex3f(system->dimension,system->dimension,system->dimension);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(0.0f,system->dimension,0.0f);
		glVertex3f(0.0f,0.0f,system->dimension);
		glVertex3f(0.0f,system->dimension,system->dimension);
		glVertex3f(system->dimension,0.0f,0.0f);
		glVertex3f(system->dimension,system->dimension,0.0f);
		glVertex3f(system->dimension,0.0f,system->dimension);
		glVertex3f(system->dimension,system->dimension,system->dimension);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(system->dimension,0.0f,0.0f);
		glVertex3f(0.0f,system->dimension,0.0f);
		glVertex3f(system->dimension,system->dimension,0.0f);
		glVertex3f(0.0f,0.0f,system->dimension);
		glVertex3f(system->dimension,0.0f,system->dimension);
		glVertex3f(0.0f,system->dimension,system->dimension);
		glVertex3f(system->dimension,system->dimension,system->dimension);
	glEnd();
}
void RenderConnectivity(SYSTEM *system)
{
	int i,j,k,l,m,n,o,p,im,in,io;
	ATOM *atom1;
	ATOM *atom2;
  //Draw Connectivity
	glBegin(GL_LINES);
		glColor3f(0.1f, 0.1f, 0.8f);
		for(i=0; i<system->nSlice; i++)
		{
			for(j=0; j<system->nSlice; j++)
			{
				for(k=0; k<system->nSlice; k++)
				{
					int idx1 = i*system->nSlice*system->nSlice+j*system->nSlice+k;
					for(l=0; l<system->gridCount[idx1]; l++)
					{
						atom1 = system->grid[idx1][l];
						for(im=i-1; im<i+2; im++)
						{
							for(in=j-1; in<j+2; in++)
							{
								for(io=k-1; io<k+2; io++)
								{
									m=im;
									n=in;
									o=io;
									if(m<0) m+=system->nSlice;
									if(n<0) n+=system->nSlice;
									if(o<0) o+=system->nSlice;
									if(m>=system->nSlice) m-=system->nSlice;
									if(n>=system->nSlice) n-=system->nSlice;
									if(o>=system->nSlice) o-=system->nSlice;
									int idx2 = m*system->nSlice*system->nSlice+n*system->nSlice+o;
									for(p=0; p<system->gridCount[idx2]; p++)
									{
										atom2 = system->grid[idx2][p];
										if(atom1->mol == atom2->mol) continue;
										if(	abs(atom1->x-atom2->x)<system->dimension/2.0 &&
												abs(atom1->y-atom2->y)<system->dimension/2.0 &&
												abs(atom1->z-atom2->z)<system->dimension/2.0 )
										{
											glVertex3f(	atom1->x,atom1->y,atom1->z);
											glVertex3f(	atom2->x,atom2->y,atom2->z);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	glEnd();
}

void RenderAtoms(SYSTEM *system)
{
	int i,j;
	ATOM *atom1;
	//Draw Atoms
	glBegin(GL_POINTS);
  	for(i=0;i<system->nAllMolecules;i++)
  	{
  		for(j=0;j<(system->allMolecules+i)->nAtoms; j++)
  		{
			atom1 = ((system->allMolecules+i)->atoms+j);	
			glColor3f(AtomColor[atom1->type][0],AtomColor[atom1->type][1],AtomColor[atom1->type][2]);
			glVertex3f(atom1->x,atom1->y,atom1->z);
		}
  	}
  glEnd();
}

void RenderBonds(SYSTEM *system)
{
	int i,j;
	ATOM *atom1;
	ATOM *atom2;

	//Draw Bonds
  glBegin(GL_LINES);
  	for(i=0;i<system->nAllMolecules;i++)
  	{
  		for(j=0;j<(system->allMolecules+i)->nBonds; j++)
  		{
  			atom1 = ((system->allMolecules+i)->bonds+j)->atom1;
			atom2 = ((system->allMolecules+i)->bonds+j)->atom2;
			if(	abs(atom1->x-atom2->x)<system->dimension/2.0 &&
  				abs(atom1->y-atom2->y)<system->dimension/2.0 &&
  				abs(atom1->z-atom2->z)<system->dimension/2.0)
  			{
					double midx = (  atom1->x + atom2->x ) / 2.0;
					double midy = (  atom1->y + atom2->y ) / 2.0;
					double midz = (  atom1->z + atom2->z ) / 2.0;
					glColor3f(AtomColor[atom1->type][0],AtomColor[atom1->type][1],AtomColor[atom1->type][2]);
					glVertex3f(	atom1->x,atom1->y,atom1->z);
					glVertex3f(	midx, midy, midz);
					glColor3f(AtomColor[atom2->type][0],AtomColor[atom2->type][1],AtomColor[atom2->type][2]);
					glVertex3f(	atom2->x,atom2->y,atom2->z);
					glVertex3f(	midx, midy, midz);
				}
				else
				{
				//FIXME	
				}
			}
  	}
  glEnd();
}
int RenderSystem(SYSTEM *system)
{
	int i=0;
	int j=0;
	int k,l,im,in,io,m,n,o,p;
	ATOM *atom1;
	ATOM *atom2;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glTranslatef(0.0f, 0.0f, -system->dimension*2.0);
    glMultMatrixf(GLWin.Transform.M);										// NEW: Apply Dynamic Transform
//  glRotatef((double)(system->step)/36.0, 1.0f, 1.0f, 0.0f);
  glTranslatef(-system->dimension/2.0, -system->dimension/2.0, -system->dimension/2.0);

	RenderCell(system);
//	RenderConnectivity(system);
	RenderAtoms(system);
	RenderBonds(system);
  
  //Draw Grid;
  glBegin(GL_LINES);
		glColor3f(0.0f, 0.7f, 0.0f);
		for(i=0; i<system->nSlice; i++)
		{
			glVertex3f(0.0f, system->dimension/system->nSlice*i, 0.0f);
			glVertex3f(system->dimension, system->dimension/system->nSlice*i, 0.0f);
		}
 		for(i=0; i<system->nSlice; i++)
		{
			glVertex3f( system->dimension/system->nSlice*i,0.0f, 0.0f);
			glVertex3f( system->dimension/system->nSlice*i,system->dimension, 0.0f);
		}
 glEnd();
  if (GLWin.doubleBuffered)
  {
      glXSwapBuffers(GLWin.dpy, GLWin.win);
  }

	return 1;
}
/* function called when our window is resized (should only happen in window mode) */
void resizeGLScene(unsigned int width, unsigned int height)
{
    if (height == 0)    /* Prevent A Divide By Zero If The Window Is Too Small */
        height = 1;
    glViewport(0, 0, width, height);    /* Reset The Current Viewport And Perspective Transformation */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//		glOrtho(-(GLfloat)width / (GLfloat)height * 10.0, (GLfloat)width / (GLfloat)height * 10.0, -10.0, 10.0,  -100.0, 100.0);
    gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 0.1f, 1000.0f);
    glMatrixMode(GL_MODELVIEW);
}

/* general OpenGL initialization function */
int initGL(void)
{
    glShadeModel(GL_SMOOTH);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    /* we use resizeGLScene once to set up our initial perspective */
    resizeGLScene(GLWin.width, GLWin.height);
    ArcBallInit(GLWin.width, GLWin.height, &(GLWin.arcball));
    Matrix3fSetIdentity(&(GLWin.LastRot));
    Matrix3fSetIdentity(&(GLWin.ThisRot));
	  Matrix4fSetIdentity(&(GLWin.Transform));		// Reset Rotation

    glFlush();
    return True;
}



/* function to release/destroy our resources and restoring the old desktop */
GLvoid killGLWindow(void)
{
    if (GLWin.ctx)
    {
        if (!glXMakeCurrent(GLWin.dpy, None, NULL))
        {
            //printf("Could not release drawing context.\n");
        }
        glXDestroyContext(GLWin.dpy, GLWin.ctx);
        GLWin.ctx = NULL;
    }
    /* switch back to original desktop resolution if we were in fs */
    if (GLWin.fs)
    {
        XF86VidModeSwitchToMode(GLWin.dpy, GLWin.screen, &GLWin.deskMode);
        XF86VidModeSetViewPort(GLWin.dpy, GLWin.screen, 0, 0);
    }
    XCloseDisplay(GLWin.dpy);
}


void UpdateTransformMatrix (long milliseconds)									// Perform Motion Updates Here
{
  if (GLWin.isRClicked)													// If Right Mouse Clicked, Reset All Rotations
  {
		Matrix3fSetIdentity(&(GLWin.LastRot));								// Reset Rotation
		Matrix3fSetIdentity(&(GLWin.ThisRot));								// Reset Rotation
	  Matrix4fSetRotationFromMatrix3f(&(GLWin.Transform), &(GLWin.ThisRot));		// Reset Rotation
  }

  if (!GLWin.isDragging)												// Not Dragging
  {
    if (GLWin.isClicked)												// First Click
    {
			GLWin.isDragging = 1;										// Prepare For Dragging
			GLWin.LastRot = GLWin.ThisRot;										// Set Last Static Rotation To Last Dynamic One
			ArcBallClick(&(GLWin.MousePt), &(GLWin.arcball));								// Update Start Vector And Prepare For Dragging
    }
  }
  else
  {
    if (GLWin.isClicked)												// Still Clicked, So Still Dragging
    {
      Quat4fT     ThisQuat;

      ArcBallDrag(&(GLWin.MousePt), &ThisQuat, &(GLWin.arcball));						// Update End Vector And Get Rotation As Quaternion
      Matrix3fSetRotationFromQuat4f(&(GLWin.ThisRot), &ThisQuat);		// Convert Quaternion Into Matrix3fT
      Matrix3fMulMatrix3f(&(GLWin.ThisRot), &(GLWin.LastRot));				// Accumulate Last Rotation Into This One
      Matrix4fSetRotationFromMatrix3f(&(GLWin.Transform), &(GLWin.ThisRot));	// Set Our Final Transform's Rotation From This One
    }
    else														// No Longer Dragging
        GLWin.isDragging = 0;
  }
}

void GraphicOutput(SYSTEM *system)
{
  XEvent event;
	struct timeval tv, tickCount;

	while(XPending(GLWin.dpy) > 0)
	{
		XNextEvent(GLWin.dpy, &event);
		switch(event.type)
		{
			case Expose:
				if(event.xexpose.count != 0)	break;
				RenderSystem(system);   	break;
			case ButtonPress:
				switch( event.xbutton.button ) 
				{
					case 1: GLWin.isClicked = 1; break;
					case 3: GLWin.isRClicked = 1; break;
				}
				break;
			case ButtonRelease:
				switch( event.xbutton.button ) 
				{
					case 1: GLWin.isClicked = 0; break;
					case 3: GLWin.isRClicked = 0; break;
				}
				break;
			case MotionNotify:
				GLWin.MousePt.s.X = event.xmotion.x;
				GLWin.MousePt.s.Y = event.xmotion.y;
				break;
			default:
				break;
		}
	}

	gettimeofday( &tv, NULL );
	tickCount.tv_sec = tv.tv_sec - GLWin.lastTickCount.tv_sec;
	tickCount.tv_usec = tv.tv_usec - GLWin.lastTickCount.tv_usec;

	UpdateTransformMatrix(tickCount.tv_usec / 1000 + tickCount.tv_sec * 1000);
	GLWin.lastTickCount = tickCount;
	RenderSystem(system);
}
/* this function creates our window and sets it up properly */
/* FIXME: bits is currently unused */
Bool createGLWindow(char* title, int width, int height, int bits,
                    Bool fullscreenflag)
{
    XVisualInfo *vi;
    Colormap cmap;
    int dpyWidth, dpyHeight;
    int i;
    int glxMajorVersion, glxMinorVersion;
    int vidModeMajorVersion, vidModeMinorVersion;
    XF86VidModeModeInfo **modes;
    int modeNum;
    int bestMode;
    Atom wmDelete;
    Window winDummy;
    unsigned int borderDummy;
    
    GLWin.fs = fullscreenflag;
    /* set best mode to current */
    bestMode = 0;
    /* get a connection */
    GLWin.dpy = XOpenDisplay(0);
    GLWin.screen = DefaultScreen(GLWin.dpy);
    XF86VidModeQueryVersion(GLWin.dpy, &vidModeMajorVersion,
        &vidModeMinorVersion);
    //printf("XF86VidModeExtension-Version %d.%d\n", vidModeMajorVersion,
    //    vidModeMinorVersion);
    XF86VidModeGetAllModeLines(GLWin.dpy, GLWin.screen, &modeNum, &modes);
    /* save desktop-resolution before switching modes */
    GLWin.deskMode = *modes[0];
    /* look for mode with requested resolution */
    for (i = 0; i < modeNum; i++)
    {
        if ((modes[i]->hdisplay == width) && (modes[i]->vdisplay == height))
        {
            bestMode = i;
        }
    }
    /* get an appropriate visual */
    vi = glXChooseVisual(GLWin.dpy, GLWin.screen, attrListDbl);
    if (vi == NULL)
    {
        vi = glXChooseVisual(GLWin.dpy, GLWin.screen, attrListSgl);
        GLWin.doubleBuffered = False;
     //   printf("Only Singlebuffered Visual!\n");
    }
    else
    {
        GLWin.doubleBuffered = True;
     //   printf("Got Doublebuffered Visual!\n");
    }
    glXQueryVersion(GLWin.dpy, &glxMajorVersion, &glxMinorVersion);
    //printf("glX-Version %d.%d\n", glxMajorVersion, glxMinorVersion);
    /* create a GLX context */
    GLWin.ctx = glXCreateContext(GLWin.dpy, vi, 0, GL_TRUE);
    /* create a color map */
    cmap = XCreateColormap(GLWin.dpy, RootWindow(GLWin.dpy, vi->screen),
        vi->visual, AllocNone);
    GLWin.attr.colormap = cmap;
    GLWin.attr.border_pixel = 0;

    if (GLWin.fs)
    {
        XF86VidModeSwitchToMode(GLWin.dpy, GLWin.screen, modes[bestMode]);
        XF86VidModeSetViewPort(GLWin.dpy, GLWin.screen, 0, 0);
        dpyWidth = modes[bestMode]->hdisplay;
        dpyHeight = modes[bestMode]->vdisplay;
        //printf("Resolution %dx%d\n", dpyWidth, dpyHeight);
        XFree(modes);
    
        /* create a fullscreen window */
        GLWin.attr.override_redirect = True;
        GLWin.attr.event_mask = ExposureMask | KeyPressMask | ButtonPressMask | 
            PointerMotionMask | ButtonReleaseMask | StructureNotifyMask ;
        GLWin.win = XCreateWindow(GLWin.dpy, RootWindow(GLWin.dpy, vi->screen),
            0, 0, dpyWidth, dpyHeight, 0, vi->depth, InputOutput, vi->visual,
            CWBorderPixel | CWColormap | CWEventMask | CWOverrideRedirect,
            &GLWin.attr);
        XWarpPointer(GLWin.dpy, None, GLWin.win, 0, 0, 0, 0, 0, 0);
		XMapRaised(GLWin.dpy, GLWin.win);
        XGrabKeyboard(GLWin.dpy, GLWin.win, True, GrabModeAsync,
            GrabModeAsync, CurrentTime);
        XGrabPointer(GLWin.dpy, GLWin.win, True, ButtonPressMask,
            GrabModeAsync, GrabModeAsync, GLWin.win, None, CurrentTime);
    }
    else
    {
        /* create a window in window mode*/
        GLWin.attr.event_mask = ExposureMask | KeyPressMask | ButtonPressMask |
            PointerMotionMask | ButtonReleaseMask | StructureNotifyMask;
        GLWin.win = XCreateWindow(GLWin.dpy, RootWindow(GLWin.dpy, vi->screen),
            0, 0, width, height, 0, vi->depth, InputOutput, vi->visual,
            CWBorderPixel | CWColormap | CWEventMask, &GLWin.attr);
        /* only set window title and handle wm_delete_events if in windowed mode */
        wmDelete = XInternAtom(GLWin.dpy, "WM_DELETE_WINDOW", True);
        XSetWMProtocols(GLWin.dpy, GLWin.win, &wmDelete, 1);
        XSetStandardProperties(GLWin.dpy, GLWin.win, title,
            title, None, NULL, 0, NULL);
        XMapRaised(GLWin.dpy, GLWin.win);
    }       
    /* connect the glx-context to the window */
    glXMakeCurrent(GLWin.dpy, GLWin.win, GLWin.ctx);
    XGetGeometry(GLWin.dpy, GLWin.win, &winDummy, &GLWin.x, &GLWin.y,
        &GLWin.width, &GLWin.height, &borderDummy, &GLWin.depth);
    //printf("Depth %d\n", GLWin.depth);
    //if (glXIsDirect(GLWin.dpy, GLWin.ctx)) 
        //printf("Congrats, you have Direct Rendering!\n");
    //else
        //printf("Sorry, no Direct Rendering possible!\n");
     initGL();
    return True;    
}


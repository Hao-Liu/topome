#include "graphic.h"
#include <stdio.h>

int main(int argc, char **argv)
{
  // default to graphic mode
  int graphic = 1;
  
  // default to verbose mode
  int verbose = 1;

	// system structure for global variables
  SYSTEM system;
  
  // X event handle
  GLWin.isClicked=0;
  GLWin.isRClicked=0;
  GLWin.isDragging=0;

  if(!InitSystem(&system, argc, argv)) return -1;
  CreateMolecules(&system);


	// relax to min energy mode
	printf("Relaxing...\n");
	for(;system.dimension>system.dimensionTarget;system.dimension*=0.999)
	{
		RelaxMolecules(&system);
	}
	ResetVelocity(&system);

  // create a X window if graphic mode is on
  if(graphic)
  {
		// default to windowed mode
		GLWin.fs = False;
		createGLWindow("TOPOME", 640, 480, 24, GLWin.fs);
  }
  
	// iterating simulation
	for(system.step=0; system.step<system.nSteps; system.step++)
	{
		UpdateMolecules(&system);
		if(verbose && system.step%system.verboseInterval==0)
		{
			VerboseOutput(&system);
		}
		if(graphic && system.step%system.graphicInterval==0)
		{
			GraphicOutput(&system);
		}
	}
	
	// cleaning up 
  ReleaseMolecules(&system);
  killGLWindow();
	return 0;
}

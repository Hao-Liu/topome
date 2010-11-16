#include "graphic.h"
#include <stdio.h>

int main(int argc, char **argv)
{
  /* default to windowed mode */
  GLWin.fs = False;
  
  
  SYSTEM system;
  
  if(argc!=2)
  {
  	printf("Invalid arguments!\n");
  	return 1;
  }
  InitSystem(&system, argv[1]);
  
  createGLWindow("My MD Test", 640, 480, 24, GLWin.fs);
	CreateMolecules(&system);

	for(;system.dimension>system.dimensionTarget;system.dimension*=0.999)
	{
		RelaxMolecules(&system);
	}
	ResetVelocity(&system);
	for(system.step=0; system.step<system.nSteps; system.step++)
	{
		UpdateMolecules(&system);
		if(system.step%100 == 0)
		{
		printf("%d %2.5lf %2.2lf %2.2lf %2.2lf\n",system.step, system.dimension, system.potentialEnergy, system.kineticEnergy, 
			system.potentialEnergy+system.kineticEnergy);
//			usleep(10000);
		}
//		RenderMolecules(&system);
	}
  ReleaseMolecules(&system);
  killGLWindow();
	return 0;
}

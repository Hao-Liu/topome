#include "graphic.h"
#include <stdio.h>

int main(int argc, char **argv)
{
  /* default to windowed mode */
  GLWin.fs = False;
  
  int nSteps=500000;
  
  SYSTEM system;
  
  if(argc!=2)
  {
  	printf("Invalid arguments!\n");
  	return 1;
  }
  InitSystem(&system, argv[1]);
  
  createGLWindow("My MD Test", 640, 480, 24, GLWin.fs);
	CreateMolecules(&system);

	int i=0;
	for(i=0; i<nSteps; i++)
	{
		system.step = i;
		UpdateMolecules(&system);
		if(i%100 == 0)
		{
		printf("%2.5lf %2.5e  %2.2e %2.2e %2.2e\n",system.r2min, system.vMax, system.potentialEnergy, system.kineticEnergy, 
			system.potentialEnergy+system.kineticEnergy);
//			usleep(10000);
		RenderMolecules(&system);
		}
	}
  ReleaseMolecules(&system);
  killGLWindow();
	return 0;
}

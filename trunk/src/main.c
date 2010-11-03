#include "graphic.h"
#include <stdio.h>

int main(int argc, char **argv)
{
  /* default to windowed mode */
  GLWin.fs = False;
  
  int nSteps=10000;
  
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
		UpdateMolecules(&system);
		if(i%100 == 0)
		printf("%7.2e %2.2lf %2.2lf %2.2lf\n", system.r2min, system.potentialEnergy, system.kineticEnergy, 
			system.potentialEnergy+system.kineticEnergy);
//			sleep(1);
		RenderMolecules(&system);
	}
  ReleaseMolecules(&system);
  killGLWindow();
	return 0;
}

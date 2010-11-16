#include "verbose.h"
void VerboseOutput(SYSTEM *system)
{
	printf("%d %2.2lf %2.2lf %2.2lf\n",system->step, system->potentialEnergy, system->kineticEnergy, 
	system->potentialEnergy+system->kineticEnergy);
}

#include <stdio.h>
#include "system.h"

void 
verbose_output (System *system)
{
	printf ("%d %2.2lf %2.2lf %2.2lf\n",system->step, 
	        system->potential_energy, system->kinetic_energy, 
	        system->potential_energy + system->kinetic_energy);
}

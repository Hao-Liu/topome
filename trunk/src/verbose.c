#include <stdio.h>
#include "system.h"

void 
verbose_output (System *system)
{
	printf ("%d %2.2lf %2.2lf %2.2lf\t",system->step, 
	        system->potential_energy, system->kinetic_energy, 
	        system->potential_energy + system->kinetic_energy);
	printf ("%d\n", system->number_atom);
}

#ifndef _tpm_force_h
#define _tpm_force_h

#include "types.h"

void 
get_minimun_image (double *dx, double *dy, double *dz, System *tpm_system);
void
calculate_bond_force (System *tpm_system);
void
calculate_angle_force (System *tpm_system);
void
calculate_pair_force_grid (System *tpm_system);
void
calculate_pair_force_all (System *tpm_system);
void
release_grid_list (System *tpm_system);
void
create_grid_list (System *tpm_system);
void
calculate_aom_force (System *tpm_system);
void
calculate_force (System *tpm_system);
void
integrate (int relax, System *tpm_system);
void
calculate_kinetic_energy (System *tpm_system);
#endif //_tpm_molecule_h

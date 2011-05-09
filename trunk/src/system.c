/*
 * /__ __X  _ \/  __\/  _ \/ \__/|/  __/ Copyright (C) 2010 - Liu Hao
 *   / \ | / \||  \/|| / \|| |\/|||  \  
 *   | | | \_/||  __/| \_/|| |  |||  /_  system.c
 *   \_/ \____/\_/   \____/\_/  \|\____\
 * 
 * global environment initiation, manipulation, and cleaning.
 *
 * This file is part of topome - a molecular dynamics and 
 * topology analysis package for coordination systems
 *
 * topome is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * topome is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with topome; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, 
 * Boston, MA  02110-1301  USA
 */


#include <errno.h>
#include <error.h>
#include <string.h>
#include <libintl.h>
#include <stdlib.h>
#include <stdio.h>
#include "force.h"
#include "system.h"
#include "graphic.h"
#include "parse.h"
#include "verbose.h"

void
init_atom (int *mol_index, MOLECULETYPE *mol_type, MOLECULE *mol, 
                    System *tpm_system)
{
  int i;
  ATOM *atom;
  ATOM *atom_type;

  //Initiate molecules at random position
  double molecular_position_x = (double)rand() / (double)RAND_MAX *
                                tpm_system->dimension;
  double molecular_position_y = (double)rand() / (double)RAND_MAX *
                                tpm_system->dimension;
  double molecular_position_z = (double)rand() / (double)RAND_MAX *
                                tpm_system->dimension;

	mol->atoms = malloc (mol->nAtoms * sizeof (ATOM));
  
  for (i=0; i<mol->nAtoms; i++)
  {
    atom = mol->atoms+i;
    atom_type = mol_type->molecule.atoms+i;
    
    atom->x = atom_type->x + molecular_position_x; 
    atom->y = atom_type->y + molecular_position_y; 
    atom->z = atom_type->z + molecular_position_z;
    atom->oldx = atom->x; 
    atom->oldy = atom->y; 
    atom->oldz = atom->z; 

    //Init zero velocity	
    atom->vx = atom_type->vx ; 
    atom->vy = atom_type->vy ; 
    atom->vz = atom_type->vz ; 
		
    atom->mol = *mol_index;
    atom->ax = atom_type->ax ; 
    atom->ay = atom_type->ay ; 
    atom->az = atom_type->az ; 
    atom->mass = atom_type->mass ; 
    atom->type = atom_type->type ; 
    atom->charge = atom_type->charge ; 
		atom->nonbond_a = atom_type->nonbond_a ;
		atom->nonbond_b = atom_type->nonbond_b ;
  }
}

void
init_bond (MOLECULETYPE *mol_type, MOLECULE *mol, 
                    System *tpm_system)
{
  int i;
  BOND *bond;
  BOND *bond_type;

  //Initialize every bonds
  mol->bonds = malloc (mol->nBonds * sizeof (BOND));
  for (i=0; i<mol->nBonds; i++)
  {
    bond = mol->bonds+i;
    bond_type = mol_type->molecule.bonds+i;
    
	  bond->type = bond_type->type ;
	  bond->idxAtom1 = bond_type->idxAtom1 ;
	  bond->idxAtom2 = bond_type->idxAtom2 ;
		bond->r0 = bond_type->r0 ;
		bond->k2 = bond_type->k2 ;
		bond->bondType = bond_type->bondType ;
	  bond->atom1 = (mol->atoms + bond_type->idxAtom1) ;
	  bond->atom2 = (mol->atoms + bond_type->idxAtom2) ;
  }
}

void
init_angle (MOLECULETYPE *mol_type, MOLECULE *mol, 
                    System *tpm_system)
{
  int i;
  ANGLE *angle;
  ANGLE *angle_type;

  mol->angles = malloc (mol->nAngles * sizeof (ANGLE));
  for(i=0; i<mol->nAngles; i++)
  {
    angle = mol->angles+i;
    angle_type = mol_type->molecule.angles+i;
    
	  angle->type = angle_type->type ;
	  angle->theta0 = angle_type->theta0 ;
	  angle->k2 = angle_type->k2 ;
	  angle->idxAtom1 = angle_type->idxAtom1 ;
	  angle->idxAtom2 = angle_type->idxAtom2 ;
	  angle->idxAtom3 = angle_type->idxAtom3 ;
	  angle->atom1 = (mol->atoms + angle_type->idxAtom1) ;
	  angle->atom2 = (mol->atoms + angle_type->idxAtom2) ;
	  angle->atom3 = (mol->atoms + angle_type->idxAtom3) ;
  }
}

void
init_dihedral (MOLECULETYPE *mol_type, MOLECULE *mol, 
                    System *tpm_system)
{
  int i;
  DIHEDRAL *dihedral;
  DIHEDRAL *dihedral_type;

  mol->dihedrals = malloc (mol->nDihedrals * sizeof (DIHEDRAL));
  for(i=0; i<mol->nDihedrals; i++)
  {
    dihedral = mol->dihedrals+i;
    dihedral_type = mol_type->molecule.dihedrals+i;
    
	  dihedral->type = dihedral_type->type ;
	  dihedral->kphi = dihedral_type->kphi ;
	  dihedral->n = dihedral_type->n ;
	  dihedral->phi0 = dihedral_type->phi0 ;
	  dihedral->idxAtom1 = dihedral_type->idxAtom1 ;
	  dihedral->idxAtom2 = dihedral_type->idxAtom2 ;
	  dihedral->idxAtom3 = dihedral_type->idxAtom3 ;
	  dihedral->idxAtom4 = dihedral_type->idxAtom4 ;
	  dihedral->atom1 = (mol->atoms + dihedral_type->idxAtom1) ;
	  dihedral->atom2 = (mol->atoms + dihedral_type->idxAtom2) ;
	  dihedral->atom3 = (mol->atoms + dihedral_type->idxAtom3) ;
	  dihedral->atom4 = (mol->atoms + dihedral_type->idxAtom4) ;
  }
}

void
init_improper (MOLECULETYPE *mol_type, MOLECULE *mol, 
                    System *tpm_system)
{
  int i;
  IMPROPER *improper;
  IMPROPER *improper_type;

  mol->impropers = malloc (mol->nImpropers * sizeof (IMPROPER));
  for (i=0; i<mol->nImpropers; i++)
  {
    improper = mol->impropers+i;
    improper_type = mol_type->molecule.impropers+i;
    
	  improper->type = improper_type->type ;
	  improper->kchi = improper_type->kchi ;
	  improper->n = improper_type->n ;
	  improper->chi0 = improper_type->chi0 ;
	  improper->idxAtom1 = improper_type->idxAtom1 ;
	  improper->idxAtom2 = improper_type->idxAtom2 ;
	  improper->idxAtom3 = improper_type->idxAtom3 ;
	  improper->idxAtom4 = improper_type->idxAtom4 ;
	  improper->atom1 = (mol->atoms + improper_type->idxAtom1) ;
	  improper->atom2 = (mol->atoms + improper_type->idxAtom2) ;
	  improper->atom3 = (mol->atoms + improper_type->idxAtom3) ;
	  improper->atom4 = (mol->atoms + improper_type->idxAtom4) ;
  }
}

void
init_molecule_type (int *mol_index, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i;

  for (i=0; i<mol_type->nMolecules; i++)
  {      
  	MOLECULE *mol = tpm_system->molecule + *(mol_index);
  	
	  mol->nAtoms = mol_type->molecule.nAtoms;
	  mol->nBonds = mol_type->molecule.nBonds;
	  mol->nAngles = mol_type->molecule.nAngles; 
	  mol->nDihedrals = mol_type->molecule.nDihedrals;
	  mol->nImpropers = mol_type->molecule.nImpropers;
			
	
    init_atom (mol_index, mol_type, mol, tpm_system);
    init_bond (mol_type, mol, tpm_system);
    init_angle (mol_type, mol, tpm_system);
    init_dihedral (mol_type, mol, tpm_system);
    init_improper (mol_type, mol, tpm_system);

	  (*mol_index)++;
  }
}

void 
create_molecule (System *tpm_system)
{
	int i=0;
	
	//get total number of molecules in the system
	tpm_system->number_molecule = 0;
	for (i=0; i<tpm_system->number_molecule_type; i++)
	{
		tpm_system->number_molecule += (tpm_system->molecule_type+i)->nMolecules;
	}
	
	//allocate memory for all molecules
	tpm_system->molecule = malloc (tpm_system->number_molecule * 
	                               sizeof (MOLECULE));
	
	//initialize molecules
	int mol_index=0;
	for (i=0; i<tpm_system->number_molecule_type; i++)
	{
	  MOLECULETYPE *mol_type = tpm_system->molecule_type+i;
	  
    init_molecule_type (&mol_index, mol_type, tpm_system);
	}
}

void 
init_system (System *tpm_system, char *input_file)
{
	FILE *fp = fopen (input_file, "r");

  if (fp == NULL)
  {
  	error (EXIT_FAILURE, errno, gettext ("Can't open input file %s"), input_file);
  }
  
  //read all information in the input file
  parse_input_file (fp, tpm_system);

  //initialize all the molecules in the system
	create_molecule (tpm_system);
}

void	
update_system (int relax, System *tpm_system)
{
	integrate (relax, tpm_system);
	create_grid_list (tpm_system);
	calculate_force (tpm_system);
	release_grid_list (tpm_system);
	calculate_kinetic_energy (tpm_system);
}



void
reset_velosity(System *tpm_system)
{
	int i,j;
	for(i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule+i;
		for(j=0; j<mol->nAtoms; j++)
		{
		  ATOM *atom = mol->atoms+j;
			atom->oldx=atom->x;
			atom->oldy=atom->y;
			atom->oldz=atom->z;
		}
	}
}


void 
relax_system (System *tpm_system)
{
  fprintf (stdout, gettext ("Relaxing...\n"));
  for (; tpm_system->dimension > tpm_system->relaxed_dimension;
       tpm_system->dimension *= tpm_system->relax_ratio)
  {
    tpm_system->half_dimension = tpm_system->dimension*0.5;
    update_system (1, tpm_system);
  }
  tpm_system->half_dimension = tpm_system->dimension*0.5;
  reset_velosity (tpm_system);
}

void
run_system (System *tpm_system, int verbose_enabled, int graphic_enabled, 
            GLWindow *gl_window)
{
  fprintf (stdout, gettext ("Running...\n"));
  for (tpm_system->step=0; tpm_system->step<tpm_system->number_step; 
       tpm_system->step++)
  {
    update_system (0, tpm_system);
    if(verbose_enabled && 
       (tpm_system->step%tpm_system->verbose_interval == 0))
    {
	    verbose_output (tpm_system);
    }
    if(graphic_enabled && 
       (tpm_system->step%tpm_system->graphic_interval == 0))
    {
	    graphic_output (tpm_system, gl_window);
    }
  }
}

void
release_system (System *tpm_system)
{
	int i=0;
	MOLECULE *mol;

	//Release contents of every molecule
	for(i=0; i<tpm_system->number_molecule; i++)
	{
	  mol = tpm_system->molecule+i;
	  
		free (mol->atoms);
		free (mol->bonds);
		free (mol->angles);
		free (mol->dihedrals);
		free (mol->impropers);
		mol->atoms     = NULL;
		mol->bonds     = NULL;
		mol->angles    = NULL;
		mol->dihedrals = NULL;
		mol->impropers = NULL;
	}
	
	//Release molecule
	free (tpm_system->molecule);
	tpm_system->molecule = NULL;
}


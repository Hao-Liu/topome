/*
 * /__ __X  _ \/  __\/  _ \/ \__/|/  __/ Copyright (C) 2010 - Liu Hao
 *   / \ | / \||  \/|| / \|| |\/|||  \  
 *   | | | \_/||  __/| \_/|| |  |||  /_  parse.c
 *   \_/ \____/\_/   \____/\_/  \|\____\
 * 
 * parse command options and input file
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


#include <string.h>
#include <libintl.h>// for gettext()
#include <stdlib.h> // for exit()
#include <stdio.h>  // for printf()
#include <getopt.h>
#include "system.h"

#define MAJOR_VERSION 0
#define MINOR_VERSION 0

#define ATOM_TYPE_MAX 119

void 
print_banner (void)
{
  fprintf (stdout, gettext("\n /__ __X  _ \\/  __\\/  _ \\/ \\__/|/  __/ \n"));
  fprintf (stdout, gettext("   / \\ | / \\||  \\/|| / \\|| |\\/|||  \\   \n"));
  fprintf (stdout, gettext("   | | | \\_/||  __/| \\_/|| |  |||  /_  \n"));
  fprintf (stdout, gettext("   \\_/ \\____/\\_/   \\____/\\_/  \\|\\____\\\n"));
  fprintf (stdout, gettext("\nCopyright (C) 2010 Liu Hao\n\n"));
}

void
print_license (void)
{
  fprintf (stdout, gettext("License GPLv3+: GNU GPL version 3 or later"));
  fprintf (stdout, gettext(" <http://gnu.org/licenses/gpl.html>.\n"));
  fprintf (stdout, gettext("This is free software: you are free to"));
  fprintf (stdout, gettext("change and redistribute it.\n"));
  fprintf (stdout, gettext("There is NO WARRANTY, to the extent"));
  fprintf (stdout, gettext("permitted by law.\n\n"));
  fprintf (stdout, gettext("Written by Liu Hao\n"));
}

void 
print_version (void)
{
  fprintf (stdout, gettext("topome %d.%d\n"), MAJOR_VERSION, MINOR_VERSION);
  print_banner ();
  print_license ();
}

void
print_usage (void)
{
  print_banner ();
  fprintf (stdout, gettext("Usage: topome [OPTION]... FILE\n"));
  fprintf (stdout, gettext("topome: molecular dynamics and topology "));
  fprintf (stdout, gettext("analysis package for coordination system\n"));
  fprintf (stdout, gettext("Read FILE as input, FILE is mandatory and"));
  fprintf (stdout, gettext("allowed only one instance.\n\n"));
  fprintf (stdout, gettext("Mandatory arguments to long options are "));
  fprintf (stdout, gettext("mandatory for short options too.\n"));
  fprintf (stdout, gettext("  -g, --graphic            run topome in "));
  fprintf (stdout, gettext("interactive graphic mode\n"));
  fprintf (stdout, gettext("  -v, --verbose            print verbose "));
  fprintf (stdout, gettext("information\n"));
  fprintf (stdout, gettext("      --help     display this help and exit\n"));
  fprintf (stdout, gettext("      --version  output version "));
  fprintf (stdout, gettext("information and exit\n\n"));
  fprintf (stdout, gettext("Report topome issues to "));
  fprintf (stdout, gettext("<http://code.google.com/p/topome/issues>\n"));
  fprintf (stdout, gettext("Home page of topome : "));
  fprintf (stdout, gettext("<http://code.google.com/p/topome/>\n"));
  fprintf (stdout, gettext("For documentation, it's not available "));
  fprintf (stdout, gettext("yet. I'm working on it. \n"));
  fprintf (stdout, gettext("Please be patient or come over for help.\n"));
  fprintf (stdout, gettext("All kinds of help are appreciated.\n"));
}

void
parse_option (char **input_file, int *graphic_enabled, 
                    int *verbose_enabled, int argc, char** argv)
{
  int c;

  while (1) 
  {
    int option_index = 0;
    static struct option long_options[] = 
    {
      {"version", 0, 0, 0},
      {"help", 0, 0, 1},
      {"graphic", 0, 0, 'g'},
      {"verbose", 0, 0, 'v'},
      {0, 0, 0, 0}
    };

    c = getopt_long(argc, argv, "gv", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) 
    {
      case 0:
        print_version ();
        exit(0);
        break;

      case 1:
        print_usage ();
        exit(0);
        break;

      case 'v':
        *verbose_enabled = 1;
        break;

      case 'g':
        *graphic_enabled = 1;
        break;

      case '?':
        break;

      default:
        printf("?? getopt returned character code 0%o ??\n", c);
        exit(1);
    }
  }

  if (argc - optind == 1) 
  {
    *input_file = malloc (strlen (argv[optind])+1);
    strcpy (*input_file, argv[optind]);
  }
  else
  {
    print_usage ();
    exit (1);
  }
}

void
parse_input_global (FILE *fp, System *tpm_system)
{
  char buffer[200];
  
	tpm_system->verbose_interval = 1000;
  tpm_system->graphic_interval = 10;
  tpm_system->time_interval = 0.5;
  tpm_system->relax_ratio      = 0.999;
  tpm_system->expansion_ratio  = 5.0;
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%d", &(tpm_system->verbose_interval));
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%d", &(tpm_system->graphic_interval));
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%lf", &(tpm_system->time_interval));
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%lf", &(tpm_system->relax_ratio));
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%lf", &(tpm_system->expansion_ratio));
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%d", &(tpm_system->number_step));
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%lf", &(tpm_system->relaxed_dimension));
  fgets (buffer, 200, fp); 
  sscanf (buffer, "%lf", &(tpm_system->radius_cut));
}

void
parse_input_bond_parameter (FILE *fp, System *tpm_system)
{
  int i = 0;
  char buffer[200];

  for (i=0; i<tpm_system->number_bond_type; i++)
  {
    fgets (buffer, 200, fp);
    sscanf (buffer, "%lf %lf", 
            &((tpm_system->bond_type+i)->r0), 
            &((tpm_system->bond_type+i)->K));
  }
}

void
parse_input_angle_parameter (FILE *fp, System *tpm_system)
{
  int i = 0;
  char buffer[200];

  for (i=0; i<tpm_system->number_angle_type; i++)
  {
    fgets (buffer, 200, fp);
    sscanf (buffer, "%lf %lf", 
            &((tpm_system->angle_type+i)->theta0), 
            &((tpm_system->angle_type+i)->K));
  }
}

void
parse_input_dihedral_parameter (FILE *fp, System *tpm_system)
{
  int i = 0;
  char buffer[200];

  for (i=0; i<tpm_system->number_dihedral_type; i++)
  {
    fgets (buffer, 200, fp);
    sscanf (buffer, "%lf %lf", 
            &((tpm_system->dihedral_type+i)->phi0), 
            &((tpm_system->dihedral_type+i)->K));
  }
}

void
parse_input_improper_parameter (FILE *fp, System *tpm_system)
{
  int i = 0;
  char buffer[200];

  for (i=0; i<tpm_system->number_improper_type; i++)
  {
    fgets (buffer, 200, fp);
    sscanf (buffer, "%lf %lf", 
            &((tpm_system->improper_type+i)->chi0), 
            &((tpm_system->improper_type+i)->K));
  }
}


void
parse_input_force_parameter (FILE *fp, System *tpm_system)
{
  char buffer[200];

  fgets (buffer, 200, fp);
  sscanf (buffer, "%d %d %d %d", 
          &(tpm_system->number_bond_type),
          &(tpm_system->number_angle_type),
          &(tpm_system->number_dihedral_type),
          &(tpm_system->number_improper_type));
  tpm_system->bond_type     = malloc (tpm_system->number_bond_type *
                                      sizeof (BONDTYPE));
  tpm_system->angle_type    = malloc (tpm_system->number_angle_type *
                                      sizeof (ANGLETYPE));
  tpm_system->dihedral_type = malloc (tpm_system->number_dihedral_type *
                                      sizeof (DIHEDRALTYPE));
  tpm_system->improper_type = malloc (tpm_system->number_improper_type *
                                      sizeof (IMPROPERTYPE));

  parse_input_bond_parameter (fp, tpm_system);
  parse_input_angle_parameter (fp, tpm_system);
  parse_input_dihedral_parameter (fp, tpm_system);
  parse_input_improper_parameter (fp, tpm_system);
}

void
parse_input_atom_detail (FILE *fp, ATOM *atom, System *tpm_system)
{
  int i;
  char buffer[200];

  const char atom_type[ATOM_TYPE_MAX][4]=
  {
    "L",
    "H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe",
    "Cs","Ba",
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra",
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
    "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",
    "Uut","Uuq","Uup","Uuh","Uus","Uuo"    
  };


  fgets (buffer, 200, fp);
  sscanf (buffer, "%lf %lf %lf %s %lf %lf", 
          &(atom->x),
          &(atom->y),
          &(atom->z),
          atom->symbol,
          &(atom->mass),
          &(atom->charge));
          
  for (i=0; i<ATOM_TYPE_MAX; i++)
  {
    if (!strcmp (atom->symbol, atom_type[i]))
    {
      atom->type = i;
      break;
    }
  }
  
  atom->vx = 0.0;
  atom->vy = 0.0;
  atom->vz = 0.0;
  atom->ax = 0.0;
  atom->ay = 0.0;
  atom->az = 0.0;
}

void
parse_input_atom_detail_auto (FILE *fp, ATOM *atom, System *tpm_system)
{
  int i;
  char buffer[200];

  const char atom_type[ATOM_TYPE_MAX][4]=
  {
    "L",
    "H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe",
    "Cs","Ba",
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra",
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
    "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",
    "Uut","Uuq","Uup","Uuh","Uus","Uuo"    
  };


  fgets (buffer, 200, fp);
  sscanf (buffer, "%lf %lf %lf %s %lf", 
          &(atom->x),
          &(atom->y),
          &(atom->z),
          atom->symbol_md,
          &(atom->charge));
  for (i=0; i<tpm_system->forcefield.num_atom_type; i++)
  {
    AtomType *at = tpm_system->forcefield.atom_type+i; 
    if (!strcmp (atom->symbol_md, at->type))
    {
      atom->mass = at->mass;
      strcpy (atom->symbol, at->element);
    }
  } 
  
  for (i=0; i<tpm_system->forcefield.num_equivalence; i++)
  {
    Equivalence *eq = tpm_system->forcefield.equivalence+i; 
    if (!strcmp (atom->symbol_md, eq->type))
    {
      strcpy (atom->symbol_nonbond, eq->nonbond);
      strcpy (atom->symbol_bond, eq->bond);
      strcpy (atom->symbol_angle, eq->angle);
      strcpy (atom->symbol_torsion, eq->torsion);
      strcpy (atom->symbol_out_of_plane, eq->out_of_plane);
    }
  } 
  
  for (i=0; i<tpm_system->forcefield.num_nonbond_lj; i++)
  {
    NonbondLJ *nb = tpm_system->forcefield.nonbond_lj+i; 
    if (!strcmp (atom->symbol_nonbond, nb->i))
    {
      atom->nonbond_a = nb->a;
      atom->nonbond_b = nb->b;
    }
  } 
  
  for (i=0; i<ATOM_TYPE_MAX; i++)
  {
    if (!strcmp (atom->symbol, atom_type[i]))
    {
      atom->type = i;
      break;
    }
  }
  
  atom->vx = 0.0;
  atom->vy = 0.0;
  atom->vz = 0.0;
  atom->ax = 0.0;
  atom->ay = 0.0;
  atom->az = 0.0;

}

void
parse_input_atom (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i;
  
  //Initiate Atoms and Calculate Center of Gravity
  double accumulated_center_of_gravity_x = 0.0;
  double accumulated_center_of_gravity_y = 0.0;
  double accumulated_center_of_gravity_z = 0.0;
  double accumulated_mass                = 0.0;
  
  for (i=0; i<mol_type->molecule.nAtoms; i++)
  {
    ATOM *atom = mol_type->molecule.atoms+i;
    
    parse_input_atom_detail (fp, atom, tpm_system);
    
    accumulated_center_of_gravity_x += atom->x * atom->mass;
    accumulated_center_of_gravity_y += atom->y * atom->mass;
    accumulated_center_of_gravity_z += atom->z * atom->mass;
    accumulated_mass                += atom->mass;
  }
  double center_of_gravity_x = accumulated_center_of_gravity_x /
                               accumulated_mass;
  double center_of_gravity_y = accumulated_center_of_gravity_y /
                               accumulated_mass;
  double center_of_gravity_z = accumulated_center_of_gravity_z /
                               accumulated_mass;

  //Adjust Origin to Center of Gravity
  for (i=0; i<mol_type->molecule.nAtoms; i++)
  {
    ATOM *atom = mol_type->molecule.atoms+i;

    atom->x -= center_of_gravity_x;
    atom->y -= center_of_gravity_y;
    atom->z -= center_of_gravity_z;
  }
}

void
parse_input_atom_auto (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i;
  
  //Initiate Atoms and Calculate Center of Gravity
  double accumulated_center_of_gravity_x = 0.0;
  double accumulated_center_of_gravity_y = 0.0;
  double accumulated_center_of_gravity_z = 0.0;
  double accumulated_mass                = 0.0;
  
  for (i=0; i<mol_type->molecule.nAtoms; i++)
  {
    ATOM *atom = mol_type->molecule.atoms+i;
    
    parse_input_atom_detail_auto (fp, atom, tpm_system);
    
    accumulated_center_of_gravity_x += atom->x * atom->mass;
    accumulated_center_of_gravity_y += atom->y * atom->mass;
    accumulated_center_of_gravity_z += atom->z * atom->mass;
    accumulated_mass                += atom->mass;
  }
  double center_of_gravity_x = accumulated_center_of_gravity_x /
                               accumulated_mass;
  double center_of_gravity_y = accumulated_center_of_gravity_y /
                               accumulated_mass;
  double center_of_gravity_z = accumulated_center_of_gravity_z /
                               accumulated_mass;

  //Adjust Origin to Center of Gravity
  for (i=0; i<mol_type->molecule.nAtoms; i++)
  {
    ATOM *atom = mol_type->molecule.atoms+i;

    atom->x -= center_of_gravity_x;
    atom->y -= center_of_gravity_y;
    atom->z -= center_of_gravity_z;
  }
}

void
parse_input_bond (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i;
  char buffer[200];
  
  //Initiate Bonds
  for (i=0; i<mol_type->molecule.nBonds; i++)
  {
    BOND *bond = mol_type->molecule.bonds + i;
    fgets (buffer, 200, fp);
    sscanf (buffer, "%d %d %d %d ", 
            &(bond->idxAtom1),
            &(bond->idxAtom2),
            &(bond->bondType),
            &(bond->type));
    (bond->idxAtom1)--;
    (bond->idxAtom2)--;
  }
}

void
parse_input_bond_auto (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i,j;
  char buffer[200];
  
  //Initiate Bonds
  for (i=0; i<mol_type->molecule.nBonds; i++)
  {
    BOND *bond = mol_type->molecule.bonds + i;
    fgets (buffer, 200, fp);
    sscanf (buffer, "%d %d %d", 
            &(bond->idxAtom1),
            &(bond->idxAtom2),
            &(bond->bondType));
    (bond->idxAtom1)--;
    (bond->idxAtom2)--;
    
    for (j=0; j<tpm_system->forcefield.num_quadratic_bond; j++)
    {
      QuadraticBond *qb = tpm_system->forcefield.quadratic_bond+j;
      bond->atom1 = mol_type->molecule.atoms + bond->idxAtom1;
      bond->atom2 = mol_type->molecule.atoms + bond->idxAtom2;
      if ( (!strcmp(bond->atom1->symbol_bond, qb->i) 
            && !strcmp(bond->atom2->symbol_bond, qb->j)) ||
           (!strcmp(bond->atom1->symbol_bond, qb->j) 
            && !strcmp(bond->atom2->symbol_bond, qb->i)))
      {
        bond->r0 = qb->r0;
        bond->k2 = qb->k2;
      }
    }
  }

}

void
parse_input_angle (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i;
  char buffer[200];

  //Initiate Angles
  for(i=0;i<mol_type->molecule.nAngles;i++)
  {
    ANGLE *angle = mol_type->molecule.angles + i;
    fgets (buffer, 200, fp);
    sscanf (buffer, "%d %d %d %d", 
            &(angle->idxAtom1),
            &(angle->idxAtom2),
            &(angle->idxAtom3),
            &(angle->type));
    (angle->idxAtom1)--;
    (angle->idxAtom2)--;
    (angle->idxAtom3)--;
  }
}

void
get_angle (MOLECULETYPE *mol_type)
{
  int i, j, k, l;
  for (i=0; i<mol_type->molecule.nBonds; i++)
  {
    BOND *bond1 = mol_type->molecule.bonds + i;
    for (j=i+1; j<mol_type->molecule.nBonds; j++)
    {
      BOND *bond2 = mol_type->molecule.bonds + j;
      int atom[4];
      atom[0] = bond1->idxAtom1;
      atom[1] = bond1->idxAtom2;
      atom[2] = bond2->idxAtom1;
      atom[3] = bond2->idxAtom2;
      for (k=0; k<=1; k++)
      {
        for (l=2; l<=3; l++)
        {
          if (atom[k]==atom[l] && atom[1-k]!=atom[5-l])
          {
            mol_type->molecule.nAngles++;
            mol_type->molecule.angles = realloc (mol_type->molecule.angles, 
                                                 mol_type->molecule.nAngles
                                                 * sizeof (ANGLE));
            ANGLE *angle = mol_type->molecule.angles + 
                           mol_type->molecule.nAngles - 1;
            angle->idxAtom2 = atom[k];
            angle->idxAtom1 = atom[1-k];
            angle->idxAtom3 = atom[5-l];
            angle->atom1    = mol_type->molecule.atoms + angle->idxAtom1;
            angle->atom2    = mol_type->molecule.atoms + angle->idxAtom2;
            angle->atom3    = mol_type->molecule.atoms + angle->idxAtom3;
          }
        }
      }
    }
  }
}

void
parse_input_angle_auto (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i,j;
  
  get_angle (mol_type);
 
  //Initiate Angles
  for (i=0; i<mol_type->molecule.nAngles; i++)
  {
    ANGLE *angle = mol_type->molecule.angles + i;
    
    for (j=0; j<tpm_system->forcefield.num_quadratic_angle; j++)
    {
      QuadraticAngle *qa = tpm_system->forcefield.quadratic_angle+j;
      angle->atom1 = mol_type->molecule.atoms + angle->idxAtom1;
      angle->atom2 = mol_type->molecule.atoms + angle->idxAtom2;
      angle->atom3 = mol_type->molecule.atoms + angle->idxAtom3;
      if ( ((!strcmp(angle->atom1->symbol_angle, qa->i) 
            && !strcmp(angle->atom3->symbol_angle, qa->k)) ||
           (!strcmp(angle->atom1->symbol_angle, qa->k) 
            && !strcmp(angle->atom3->symbol_angle, qa->i))) &&
           !strcmp(angle->atom2->symbol_angle, qa->j))
      {
        angle->theta0 = qa->theta0;
        angle->k2 = qa->k2;
      }
    }
  }

}

void
parse_input_dihedral (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i;
  char buffer[200];

  //Initiate Dihedrals
  for(i=0;i<mol_type->molecule.nDihedrals;i++)
  {
    DIHEDRAL *dihedral = mol_type->molecule.dihedrals + i;
    fgets (buffer, 200, fp);
    sscanf (buffer, "%d %d %d %d %d", 
            &(dihedral->idxAtom1),
            &(dihedral->idxAtom2),
            &(dihedral->idxAtom3),
            &(dihedral->idxAtom4),
            &(dihedral->type));
    (dihedral->idxAtom1)--;
    (dihedral->idxAtom2)--;
    (dihedral->idxAtom3)--;
    (dihedral->idxAtom4)--;
  }
}

int
is_in_improper (MOLECULETYPE *mol_type, int a, int b, int c, int d)
{
  int i,j,k;
  int sorted[4];
  sorted[0]=a;
  sorted[1]=c;
  sorted[2]=d;
  sorted[3]=b;  //central atom
  for (i=0; i<2; i++)
  {
    for (j=i+1; j<3; j++)
    {
      if (sorted[i] > sorted[j])
      {
        int t = sorted[i];
        sorted[i] = sorted[j];
        sorted[j] = t;
      }
    }
  } 

  for (i=0; i<mol_type->molecule.nImpropers; i++)
  {
    IMPROPER *improper = mol_type->molecule.impropers + i;

    int sorted_ref[4];
    sorted_ref[0] = improper->idxAtom1;
    sorted_ref[1] = improper->idxAtom3;
    sorted_ref[2] = improper->idxAtom4;
    sorted_ref[3] = improper->idxAtom2; //central atom

    for (j=0; j<2; j++)
    {
      for (k=j+1; k<3; k++)
      {
        if (sorted_ref[j] > sorted_ref[k])
        {
          int t = sorted_ref[j];
          sorted_ref[j] = sorted_ref[k];
          sorted_ref[k] = t;
        }
      }
    }

    if (sorted[0] == sorted_ref[0] && 
        sorted[1] == sorted_ref[1] && 
        sorted[2] == sorted_ref[2] && 
        sorted[3] == sorted_ref[3] )
    {
      return 1;
    }
  }
  return 0;
}

//Check if requested improper angle exist in forcefield
int
is_improper_in_forcefield ()  //TODO
{
	return 1;
}

int
is_in_dihedral (MOLECULETYPE *mol_type, int a, int b, int c, int d)
{
  int i;
  int sorted[4];
  sorted[0]=a;
  sorted[1]=b;
  sorted[2]=c;
  sorted[3]=d;

	if ( sorted[0] > sorted[3])
	{
		int t = sorted[0];
		sorted[0] = sorted[3];
		sorted[3] = t;
		
		t = sorted[1];
		sorted[1] = sorted[2];
		sorted[2] = t;
	}
	else if ( sorted[0] == sorted[3] )
	{
		if ( sorted[1] > sorted[2])
		{
			int t = sorted[0];
			sorted[0] = sorted[3];
			sorted[3] = t;
			
			t = sorted[1];
			sorted[1] = sorted[2];
			sorted[2] = t;
		}
	}

  for (i=0; i<mol_type->molecule.nDihedrals; i++)
  {
    DIHEDRAL *dihedral = mol_type->molecule.dihedrals + i;

    int sorted_ref[4];
    sorted_ref[0] = dihedral->idxAtom1;
    sorted_ref[1] = dihedral->idxAtom2;
    sorted_ref[2] = dihedral->idxAtom3;
    sorted_ref[3] = dihedral->idxAtom4;

		if ( sorted_ref[0] > sorted_ref[3])
		{
			int t = sorted_ref[0];
			sorted_ref[0] = sorted_ref[3];
			sorted_ref[3] = t;
			
			t = sorted_ref[1];
			sorted_ref[1] = sorted_ref[2];
			sorted_ref[2] = t;
		}
		else if ( sorted_ref[0] == sorted_ref[3] )
		{
			if ( sorted_ref[1] > sorted_ref[2])
			{
				int t = sorted_ref[0];
				sorted_ref[0] = sorted_ref[3];
				sorted_ref[3] = t;
				
				t = sorted_ref[1];
				sorted_ref[1] = sorted_ref[2];
				sorted_ref[2] = t;
			}
		}

    if (sorted[0] == sorted_ref[0] && 
        sorted[1] == sorted_ref[1] && 
        sorted[2] == sorted_ref[2] && 
        sorted[3] == sorted_ref[3] )
    {
      return 1;
    }
  }
  return 0;
}

//Check if requested dihedral angle exist in forcefield
int
is_dihedral_in_forcefield ()  //TODO
{
	return 1;
}

void
get_dihedral_and_improper (MOLECULETYPE *mol_type)
{
  int i, j, k, l;
  for (i=0; i<mol_type->molecule.nAngles; i++)
  {
    ANGLE *angle1 = mol_type->molecule.angles + i;
    for (j=i+1; j<mol_type->molecule.nAngles; j++)
    {
      ANGLE *angle2 = mol_type->molecule.angles + j;
      int atom[2][3];
      atom[0][0] = angle1->idxAtom1;
      atom[0][1] = angle1->idxAtom2;
      atom[0][2] = angle1->idxAtom3;
      atom[1][0] = angle2->idxAtom1;
      atom[1][1] = angle2->idxAtom2;
      atom[1][2] = angle2->idxAtom3;
      if (atom[0][1] == atom[1][1]) // share apex
      {
        for (k=0; k<3; k+=2)
        {
          for (l=0; l<3; l+=2)      
          {
            if (atom[0][k]==atom[1][l] && atom[0][2-k]!=atom[1][2-l])  // share edge : improper
            {
              int a = 2 - k;  // unshared atom in angle1
              int b = 2 - l;  // unshared atom in angle2
              if (!is_in_improper (mol_type, atom[0][k], atom[0][1], atom[0][a], atom[1][b]))
              {
                if (is_improper_in_forcefield())
                {
                  mol_type->molecule.nImpropers++;
                  mol_type->molecule.impropers = realloc (mol_type->molecule.impropers, 
                                                       mol_type->molecule.nImpropers
                                                       * sizeof (IMPROPER));
                  IMPROPER *improper = mol_type->molecule.impropers + 
                                 mol_type->molecule.nImpropers - 1;
                  improper->idxAtom2 = atom[0][1];
                  improper->idxAtom1 = atom[0][k];
                  improper->idxAtom3 = atom[0][a];
                  improper->idxAtom4 = atom[1][b];
                }
              }
            }
          }
        }
      }
      else
      {
				for ( k=0; k<=2; k+=2)
				{
					for ( l=0; l<=2; l+=2)
					{
						// share an edge: dihedral
						if ( (atom[0][1] == atom[1][k]) && (atom[1][1] == atom[0][l]) ) 
						{
							int a = 2 - l; // unshared atom in angle1
							int b = 2 - k; // unshared atom in angle2
							
							if (!is_in_dihedral (mol_type, atom[0][a], atom[0][1], 
																						 atom[1][1], atom[1][b]))
							{
								if (is_dihedral_in_forcefield())
								{
									mol_type->molecule.nDihedrals++;
									mol_type->molecule.dihedrals = realloc (mol_type->molecule.dihedrals, 
																											 mol_type->molecule.nDihedrals
																											 * sizeof (DIHEDRAL));
									DIHEDRAL *dihedral = mol_type->molecule.dihedrals + 
																 mol_type->molecule.nDihedrals - 1;
									dihedral->idxAtom1 = atom[0][a];
									dihedral->idxAtom2 = atom[0][1];
									dihedral->idxAtom3 = atom[1][1];
									dihedral->idxAtom4 = atom[1][b];
								}
							}
						}
					}
				}
      }
    }
  }
}

//compare non-apex atom symbols of two impropers 
int 
compare_improper (char *a, char* b, char* c, char* x, char* y, char* z)
{
  int i,j;
  char* sorted_1[4];
  sorted_1[0]=a;
  sorted_1[1]=b;
  sorted_1[2]=c;
  for (i=0; i<3; i++)
  {
    for (j=i; j<3; j++)
    {
      if ( strcmp( sorted_1[i], sorted_1[j]) > 0)
      {
        char* t = sorted_1[i];
        sorted_1[i] = sorted_1[j];
        sorted_1[j] = t;
      }
    }
  }

  char* sorted_2[4];
  sorted_2[0]=a;
  sorted_2[1]=b;
  sorted_2[2]=c;
  for (i=0; i<3; i++)
  {
    for (j=i; j<3; j++)
    {
      if ( strcmp( sorted_2[i], sorted_2[j]) > 0)
      {
        char* t = sorted_2[i];
        sorted_2[i] = sorted_2[j];
        sorted_2[j] = t;
      }
    }
  }
	
	if (  !strcmp(sorted_1[0], sorted_2[0]) &&
				!strcmp(sorted_1[1], sorted_2[1]) &&
				!strcmp(sorted_1[2], sorted_2[2]))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void
parse_input_dihedral_and_improper_auto (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i,j;
  
  get_dihedral_and_improper (mol_type);
 
  //Initiate Dihedrals
  for (i=0; i<mol_type->molecule.nDihedrals; i++)
  {
    DIHEDRAL *dihedral = mol_type->molecule.dihedrals + i;
    
    for (j=0; j<tpm_system->forcefield.num_torsion_1; j++)
    {
      Torsion1 *t1 = tpm_system->forcefield.torsion_1+j;
      dihedral->atom1 = mol_type->molecule.atoms + dihedral->idxAtom1;
      dihedral->atom2 = mol_type->molecule.atoms + dihedral->idxAtom2;
      dihedral->atom3 = mol_type->molecule.atoms + dihedral->idxAtom3;
      dihedral->atom4 = mol_type->molecule.atoms + dihedral->idxAtom4;

			if (((!strcmp(dihedral->atom1->symbol_torsion, t1->i) || !strcmp(t1->i, "*")) &&
			     (!strcmp(dihedral->atom2->symbol_torsion, t1->j) || !strcmp(t1->j, "*")) &&
			     (!strcmp(dihedral->atom3->symbol_torsion, t1->k) || !strcmp(t1->k, "*")) &&
			     (!strcmp(dihedral->atom4->symbol_torsion, t1->l) || !strcmp(t1->l, "*")))||
			    ((!strcmp(dihedral->atom4->symbol_torsion, t1->i) || !strcmp(t1->i, "*")) &&
			     (!strcmp(dihedral->atom3->symbol_torsion, t1->j) || !strcmp(t1->j, "*")) &&
			     (!strcmp(dihedral->atom2->symbol_torsion, t1->k) || !strcmp(t1->k, "*")) &&
			     (!strcmp(dihedral->atom1->symbol_torsion, t1->l) || !strcmp(t1->l, "*"))))
      {
        dihedral->kphi = t1->kphi;
        dihedral->n = t1->n;
        dihedral->phi0 = t1->phi0;
      }
    }
  }

  //Initiate Impropers
  for (i=0; i<mol_type->molecule.nImpropers; i++)
  {
    IMPROPER *improper = mol_type->molecule.impropers + i;
    
    for (j=0; j<tpm_system->forcefield.num_out_of_plane; j++)
    {
      OutOfPlane *o1 = tpm_system->forcefield.out_of_plane+j;
      improper->atom1 = mol_type->molecule.atoms + improper->idxAtom1;
      improper->atom2 = mol_type->molecule.atoms + improper->idxAtom2;
      improper->atom3 = mol_type->molecule.atoms + improper->idxAtom3;
      improper->atom4 = mol_type->molecule.atoms + improper->idxAtom4;
			if (!strcmp(improper->atom2->symbol_out_of_plane, o1->j))
			{
				if (compare_improper ( improper->atom1->symbol_out_of_plane,
															 improper->atom3->symbol_out_of_plane,
															 improper->atom4->symbol_out_of_plane,
															 o1->i, o1->k, o1->l))
				{
					improper->kchi = o1->kchi;
					improper->n = o1->n;
					improper->chi0 = o1->chi0;
				}
			}
    }
  }

}
void
parse_input_improper (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  int i;
  char buffer[200];

  //Initiate Impropers
  for(i=0;i<mol_type->molecule.nImpropers;i++)
  {
    IMPROPER *improper = mol_type->molecule.impropers + i;
    fgets (buffer, 200, fp);
    sscanf (buffer, "%d %d %d %d %d", 
            &(improper->idxAtom1),
            &(improper->idxAtom2),
            &(improper->idxAtom3),
            &(improper->idxAtom4),
            &(improper->type));
    (improper->idxAtom1)--;
    (improper->idxAtom2)--;
    (improper->idxAtom3)--;
    (improper->idxAtom4)--;
  }
}


void
parse_input_molecule_type (FILE *fp, MOLECULETYPE *mol_type, System *tpm_system)
{
  char buffer[200];
  
  fgets (buffer, 200, fp);
  sscanf (buffer, "%d", &(mol_type->nMolecules));
  
  fgets (buffer, 200, fp);
  sscanf (buffer, "%d %d", 
          &(mol_type->molecule.nAtoms),
          &(mol_type->molecule.nBonds));
	mol_type->molecule.nAngles = 0 ;
	mol_type->molecule.nDihedrals = 0 ;
	mol_type->molecule.nImpropers = 0 ;

  mol_type->molecule.atoms     = malloc (mol_type->molecule.nAtoms *
                                         sizeof (ATOM));
  mol_type->molecule.bonds     = malloc (mol_type->molecule.nBonds *
                                         sizeof (BOND));
  mol_type->molecule.angles     = malloc (mol_type->molecule.nAngles *
                                         sizeof (ANGLE));
  mol_type->molecule.dihedrals     = malloc (mol_type->molecule.nDihedrals *
                                         sizeof (DIHEDRAL));
  mol_type->molecule.impropers     = malloc (mol_type->molecule.nImpropers *
                                         sizeof (IMPROPER));
/* Mannully input force field
  parse_input_atom (fp, mol_type, tpm_system);
  parse_input_bond (fp, mol_type, tpm_system);
  parse_input_angle (fp, mol_type, tpm_system);
  parse_input_dihedral (fp, mol_type, tpm_system);
  parse_input_improper (fp, mol_type, tpm_system);
*/
// Automatic grenerate force field
  parse_input_atom_auto (fp, mol_type, tpm_system);
  parse_input_bond_auto (fp, mol_type, tpm_system);
  parse_input_angle_auto (fp, mol_type, tpm_system);
  parse_input_dihedral_and_improper_auto (fp, mol_type, tpm_system);
}

void
parse_input_molecule (FILE *fp, System *tpm_system)
{
  int i;
  char buffer[200];

  fgets (buffer, 200, fp);
  sscanf (buffer, "%d", &(tpm_system->number_molecule_type));
  tpm_system->molecule_type = malloc (tpm_system->number_molecule_type *
                                      sizeof(MOLECULETYPE));
                                      
  for(i=0; i<tpm_system->number_molecule_type; i++)
  {
    MOLECULETYPE * mol_type = tpm_system->molecule_type+i;
    
    parse_input_molecule_type (fp, mol_type, tpm_system);
  }
}


void 
parse_input_file (FILE *fp, System *tpm_system)
{
  parse_input_global (fp, tpm_system);

  tpm_system->dimension = tpm_system->relaxed_dimension
                          * tpm_system->expansion_ratio;
  tpm_system->half_dimension = tpm_system->dimension*0.5;
  
  parse_input_force_parameter (fp, tpm_system);
  
  parse_input_molecule (fp, tpm_system);
}

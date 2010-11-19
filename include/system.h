#ifndef _tpm_system_h
#define _tpm_system_h

typedef struct
{
	double oldx;
	double oldy;
	double oldz;
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double ax;
	double ay;
	double az;
	double mass;
	int    mol;
	int    type;
	char   symbol[4];
	double charge;
} ATOM;
typedef struct
{
	ATOM *atom1;
	ATOM *atom2;
	int idxAtom1;
	int idxAtom2;
	int type;
	int bondType;
} BOND;
typedef struct
{
	ATOM *atom1;
	ATOM *atom2;
	ATOM *atom3;
	int idxAtom1;
	int idxAtom2;
	int idxAtom3;
	int type;
} ANGLE;
typedef struct
{
	ATOM *atom1;
	ATOM *atom2;
	ATOM *atom3;
	ATOM *atom4;
	int idxAtom1;
	int idxAtom2;
	int idxAtom3;
	int idxAtom4;
	int type;
} DIHEDRAL;
typedef struct
{
	ATOM *atom1;
	ATOM *atom2;
	ATOM *atom3;
	ATOM *atom4;
	int idxAtom1;
	int idxAtom2;
	int idxAtom3;
	int idxAtom4;
	int type;
} IMPROPER;

typedef struct 
{
	int nAtoms;
	int nBonds;
	int nAngles;
	int nDihedrals;
	int nImpropers;
	ATOM *atoms;
	BOND *bonds;
	ANGLE *angles;
	DIHEDRAL *dihedrals;
	IMPROPER *impropers;
} MOLECULE;

typedef struct
{
	double K;
	double r0;
} BONDTYPE;
typedef struct
{
	double K;
	double theta0;
} ANGLETYPE;
typedef struct
{
	double K;
	double phi0;
} DIHEDRALTYPE;
typedef struct
{
	double K;
	double chi0;
} IMPROPERTYPE;

typedef struct
{
	MOLECULE molecule;
	int nMolecules;
} MOLECULETYPE;

typedef struct
{
	int step;
	int number_step;
	int verbose_interval;
	int graphic_interval;
	double dimension;
	double half_dimension;
	double relaxed_dimension;
	double relax_ratio;
	double expansion_ratio;
	int number_slice;
	int number_grid;
	double radius_cut;
	ATOM ***grid;
	int *number_atom_in_grid;
	double potential_energy;
	double kinetic_energy;

	int number_bond_type;
	BONDTYPE *bond_type;
	int number_angle_type;
	ANGLETYPE *angle_type;
	int number_dihedral_type;
	DIHEDRALTYPE *dihedral_type;
	int number_improper_type;
	IMPROPERTYPE *improper_type;
	
	int number_molecule_type;
	MOLECULETYPE *molecule_type;
	int number_molecule;
	MOLECULE *molecule;
} System;

void 
init_system (System *tpm_system, char *input_file);
void 
relax_system(System *tpm_system);
void
release_system (System *tpm_system);
//void
//run_system(System *tpm_system, int verbose_enabled, int graphic_enabled, GLWindow *gl_window);

#endif //_tpm_system_h

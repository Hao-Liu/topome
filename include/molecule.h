#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
	int    type;
	double charge;
} ATOM;
typedef struct
{
	ATOM *atom1;
	ATOM *atom2;
	int idxAtom1;
	int idxAtom2;
	int type;
} BOND;
typedef struct
{
	ATOM *atom1;
	ATOM *atom2;
	ATOM *atom3;
	int idxAtom1;
	int idxAtom2;
	int idxAtom3;
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
} BONDTYPE;
typedef struct
{
} ANGLETYPE;
typedef struct
{
} DIHEDRALTYPE;
typedef struct
{
} IMPROPERTYPE;

typedef struct
{
	MOLECULE molecule;
	int nMolecules;
} MOLECULETYPE;

typedef struct
{
	double dimension;
	double potentialEnergy;
	double kineticEnergy;
	double r2min; //DEBUG FIXME

	int nBondTypes;
	BONDTYPE bondTypes;
	int nAngleTypes;
	ANGLETYPE angleTypes;
	int nDihedralTypes;
	DIHEDRALTYPE dihedralTypes;
	int nImproperTypes;
	IMPROPERTYPE improperTypes;
	
	int nMoleculeTypes;
	MOLECULETYPE *moleculeTypes;
	int nAllMolecules;
	MOLECULE *allMolecules;
} SYSTEM;
int InitSystem(SYSTEM *system, char *infile);
int CreateMolecules(SYSTEM *system);
int	UpdateMolecules(SYSTEM *system);
int ReleaseMolecules(SYSTEM *system);

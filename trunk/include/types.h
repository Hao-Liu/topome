#ifndef _tpm_types_h
#define _tpm_types_h

typedef struct ForceFunction
{
  char version[10];
  int  reference;
  char function[50];
  char label1[20];
  char label2[20];
}ForceFunction;

typedef struct AtomType
{
  char version[10];
  int  reference;
  char type[10];
  float mass;
  char element[5];
  int  connection;
}AtomType;

typedef struct ForceType
{
  char name[200];
  int  num_function;
  ForceFunction *function;
}ForceType;

typedef struct Equivalence
{
  char version[10];
  int  reference;
  char type[10];
  char nonbond[10];
  char bond[10];
  char angle[10];
  char torsion[10];
  char out_of_plane[10];
}Equivalence;

typedef struct AutoEquivalence
{
  char version[10];
  int  reference;
  char type[10];
  char nonbond[10];
  char bond_increment[10];
  char bond[10];
  char angle_end_atom[10];
  char angle_apex_atom[10];
  char torsion_end_atom[10];
  char torsion_center_atom[10];
  char out_of_plane_end_atom[10];
  char out_of_plane_center_atom[10];
}AutoEquivalence;

typedef struct HBondDefinition
{
  char version[10];
  int  reference;
  float distance;
  float angle;
  char donor[4][10];
  char acceptor[5][10];
}HBondDefinition;

typedef struct MorseBond
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  float r0;
  float d;
  float alpha;
}MorseBond;

typedef struct QuadraticBond
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  float r0;
  float k2;
}QuadraticBond;

typedef struct QuadraticAngle
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  float theta0;
  float k2;
}QuadraticAngle;

typedef struct BondBond
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  float kbb;
}BondBond;

typedef struct BondAngle
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  float kbt1;
  float kbt2;
}BondAngle;

typedef struct Torsion1
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  char l[10];
  float kphi;
  float n;
  float phi0;
}Torsion1;

typedef struct AngleAngleTorsion1
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  char l[10];
  float kaat;
}AngleAngleTorsion1;

typedef struct OutOfPlane
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  char l[10];
  float kchi;
  float n;
  float chi0;
}OutOfPlane;

typedef struct OutOfPlaneOutOfPlane
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  char l[10];
  float koo;
}OutOfPlaneOutOfPlane;

typedef struct AngleAngle
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  char k[10];
  char l[10];
  float kaa;
}AngleAngle;

typedef struct NonbondLJ
{
  char version[10];
  int  reference;
  char i[10];
  float a;
  float b;
}NonbondLJ;

typedef struct BondIncrements
{
  char version[10];
  int  reference;
  char i[10];
  char j[10];
  float dij;
  float dji;
}BondIncrements;

typedef struct ForceField
{
  char name[200];
  int  num_force_type;
  ForceType *force_type;
  int  num_atom_type;
  AtomType *atom_type;
  int  num_equivalence;
  Equivalence *equivalence;
  int  num_auto_equivalence;
  AutoEquivalence *auto_equivalence;
  HBondDefinition *hbond_definition;
  int  num_morse_bond;
  MorseBond *morse_bond;
  int  num_quadratic_bond;
  QuadraticBond *quadratic_bond;
  int  num_quadratic_angle;
  QuadraticAngle *quadratic_angle;
  int  num_bond_bond;
  BondBond *bond_bond;
  int  num_bond_angle;
  BondAngle *bond_angle;
  int  num_torsion_1;
  Torsion1 *torsion_1;
  int  num_angle_angle_torsion_1;
  AngleAngleTorsion1 *angle_angle_torsion_1;
  int  num_out_of_plane;
  OutOfPlane *out_of_plane;
  int  num_out_of_plane_out_of_plane;
  OutOfPlaneOutOfPlane *out_of_plane_out_of_plane;
  int  num_angle_angle;
  AngleAngle *angle_angle;
  int  num_nonbond_lj;
  NonbondLJ *nonbond_lj;
  int  num_bond_increments;
  BondIncrements *bond_increments;
}ForceField;

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
	char   symbol_md[6];
	char   symbol_nonbond[6];
	char   symbol_bond[6];
	char   symbol_angle[6];
	char   symbol_torsion[6];
	char   symbol_out_of_plane[6];
	double charge;
  double nonbond_a;
  double nonbond_b;
	int    global_idx;
} ATOM;

typedef struct
{
	ATOM *atom1;
	ATOM *atom2;
	int idxAtom1;
	int idxAtom2;
	int type;
  double r0;
  double k2;
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
  double theta0;
  double k2;
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
  double kphi;
  double n;
  double phi0;
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
	double kchi;
	double n;
	double chi0;
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
	double time_interval;
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
	int number_atom;
	ATOM **atom;
  ForceField forcefield;
} System;

#include <GL/glx.h>
#include <X11/extensions/xf86vmode.h>
#include "arcball.h"
/* stuff about our window grouped together */
typedef struct {
    Display *dpy;
    int screen;
    Window win;
    GLXContext ctx;
    XSetWindowAttributes attr;
    Bool fs;
    Bool doubleBuffered;
    XF86VidModeModeInfo deskMode;
    int x, y;
    
    int isClicked, isRClicked, isDragging;
    Matrix3fT   LastRot, ThisRot;
    Matrix4fT   Transform;
    Point2fT    MousePt;
    ArcBall			arcball;
    
    unsigned int width, height;
    unsigned int depth;    
		struct timeval lastTickCount;				// Tick Counter
} GLWindow;
#endif //_tpm_types_h

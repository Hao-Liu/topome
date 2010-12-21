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



#include "molecule.h"
#define ATOM_TYPE_MAX 21

int InitSystem(SYSTEM *system, int argc, char **argv)
{
	int i=0;
	int j=0;
	int k=0;
	
	if(argc!=2)
  {
  	printf("Invalid arguments!\n");
  	return 0;
  }
  

	const char AtomType[ATOM_TYPE_MAX][4]=
	{
		"X",
		"H","He","Li","Be","B","C","N","O","F","Ne",
		"Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
	};
	FILE *fp = fopen(argv[1],"r");
	char buffer[200];
	char type[4];
	int idummy;
	system->verboseInterval = 100;
	system->graphicInterval = 10;

	fgets(buffer, 200, fp);
	sscanf(buffer, "%d %lf %lf", &(system->nSteps), &(system->dimensionTarget), &(system->rCut));

	system->dimension = system->dimensionTarget*5.0;
	
	fgets(buffer, 200, fp);
		sscanf(buffer, "%d %d %d %d", 
								&(system->nBondTypes),
								&(system->nAngleTypes),
								&(system->nDihedralTypes),
								&(system->nImproperTypes));
	system->bondTypes = malloc(system->nBondTypes*sizeof(BONDTYPE));
	system->angleTypes = malloc(system->nAngleTypes*sizeof(ANGLETYPE));
	system->dihedralTypes = malloc(system->nDihedralTypes*sizeof(DIHEDRALTYPE));
	system->improperTypes = malloc(system->nImproperTypes*sizeof(IMPROPERTYPE));
	for(i=0; i<system->nBondTypes; i++)
	{
		fgets(buffer, 200, fp);
		sscanf(buffer, "%lf %lf", 
			&((system->bondTypes+i)->r0), 
			&((system->bondTypes+i)->K));
	}
	for(i=0; i<system->nAngleTypes; i++)
	{
		fgets(buffer, 200, fp);
		sscanf(buffer, "%lf %lf", 
			&((system->angleTypes+i)->theta0), 
			&((system->angleTypes+i)->K));
	}
	for(i=0; i<system->nDihedralTypes; i++)
	{
		fgets(buffer, 200, fp);
		sscanf(buffer, "%lf %lf", 
			&((system->dihedralTypes+i)->phi0), 
			&((system->dihedralTypes+i)->K));
	}
	for(i=0; i<system->nImproperTypes; i++)
	{
		fgets(buffer, 200, fp);
		sscanf(buffer, "%lf %lf", 
			&((system->improperTypes+i)->chi0), 
			&((system->improperTypes+i)->K));
	}
	
	fgets(buffer, 200, fp);
	sscanf(buffer, "%d", &(system->nMoleculeTypes));
	system->moleculeTypes = malloc(system->nMoleculeTypes*sizeof(MOLECULETYPE));
	for(i=0; i<system->nMoleculeTypes; i++)
	{
		fgets(buffer, 200, fp);
		sscanf(buffer, "%d", &((system->moleculeTypes+i)->nMolecules));
		fgets(buffer, 200, fp);
		sscanf(buffer, "%d %d %d %d %d", 
								&((system->moleculeTypes+i)->molecule.nAtoms),
								&((system->moleculeTypes+i)->molecule.nBonds),
								&((system->moleculeTypes+i)->molecule.nAngles),
								&((system->moleculeTypes+i)->molecule.nDihedrals),
								&((system->moleculeTypes+i)->molecule.nImpropers));
		(system->moleculeTypes+i)->molecule.atoms = 
				malloc((system->moleculeTypes+i)->molecule.nAtoms*sizeof(ATOM));
		(system->moleculeTypes+i)->molecule.bonds = 
				malloc((system->moleculeTypes+i)->molecule.nBonds*sizeof(BOND));
		(system->moleculeTypes+i)->molecule.angles = 
				malloc((system->moleculeTypes+i)->molecule.nAngles*sizeof(ANGLE));
		(system->moleculeTypes+i)->molecule.dihedrals = 
				malloc((system->moleculeTypes+i)->molecule.nDihedrals*sizeof(DIHEDRAL));
		(system->moleculeTypes+i)->molecule.impropers = 
				malloc((system->moleculeTypes+i)->molecule.nImpropers*sizeof(IMPROPER));

		//Initiate Atoms and Calculate Center of Gravity
		double accCOGx=0.0;
		double accCOGy=0.0;
		double accCOGz=0.0;
		double accMass=0.0;
		for(j=0;j<(system->moleculeTypes+i)->molecule.nAtoms;j++)
		{
			fgets(buffer, 200, fp);
			sscanf(buffer, "%lf %lf %lf %s %lf %lf", 
								&(((system->moleculeTypes+i)->molecule.atoms+j)->x),
								&(((system->moleculeTypes+i)->molecule.atoms+j)->y),
								&(((system->moleculeTypes+i)->molecule.atoms+j)->z),
								type,
								&(((system->moleculeTypes+i)->molecule.atoms+j)->mass),
								&(((system->moleculeTypes+i)->molecule.atoms+j)->charge));
			for(k=0; k<ATOM_TYPE_MAX; k++)
			{
				if(!strcmp(type, AtomType[k]))
				{
					((system->moleculeTypes+i)->molecule.atoms+j)->type = k;
					break;
				}
			}
			((system->moleculeTypes+i)->molecule.atoms+j)->vx = 0.0;
			((system->moleculeTypes+i)->molecule.atoms+j)->vy = 0.0;
			((system->moleculeTypes+i)->molecule.atoms+j)->vz = 0.0;
			((system->moleculeTypes+i)->molecule.atoms+j)->ax = 0.0;
			((system->moleculeTypes+i)->molecule.atoms+j)->ay = 0.0;
			((system->moleculeTypes+i)->molecule.atoms+j)->az = 0.0;
			
			accCOGx+=(((system->moleculeTypes+i)->molecule.atoms+j)->x)*
								(((system->moleculeTypes+i)->molecule.atoms+j)->mass);
			accCOGy+=(((system->moleculeTypes+i)->molecule.atoms+j)->y)*
								(((system->moleculeTypes+i)->molecule.atoms+j)->mass);
			accCOGz+=(((system->moleculeTypes+i)->molecule.atoms+j)->z)*
								(((system->moleculeTypes+i)->molecule.atoms+j)->mass);
			accMass+=(((system->moleculeTypes+i)->molecule.atoms+j)->mass);
		}
		double COGx=accCOGx/accMass;
		double COGy=accCOGy/accMass;
		double COGz=accCOGz/accMass;
		
		//Adjust Origin to Center of Gravity
		for(j=0;j<(system->moleculeTypes+i)->molecule.nAtoms;j++)
		{
			(((system->moleculeTypes+i)->molecule.atoms+j)->x)-=COGx;
			(((system->moleculeTypes+i)->molecule.atoms+j)->y)-=COGy;
			(((system->moleculeTypes+i)->molecule.atoms+j)->z)-=COGz;
		}
		
		//Initiate Bonds
		for(j=0;j<(system->moleculeTypes+i)->molecule.nBonds;j++)
		{
			fgets(buffer, 200, fp);
			sscanf(buffer, "%d %d %d %d ", 
								&(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom1),
								&(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom2),
								&(((system->moleculeTypes+i)->molecule.bonds+j)->bondType),
								&(((system->moleculeTypes+i)->molecule.bonds+j)->type));
			(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom1)--;
			(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom2)--;
		}

		//Initiate Angles
		for(j=0;j<(system->moleculeTypes+i)->molecule.nAngles;j++)
		{
			fgets(buffer, 200, fp);
			sscanf(buffer, "%d %d %d %d", 
								&(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom1),
								&(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom2),
								&(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom3),
								&(((system->moleculeTypes+i)->molecule.angles+j)->type));
			(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom1)--;
			(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom2)--;
			(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom3)--;
		}

		//Initiate Dihedrals
		for(j=0;j<(system->moleculeTypes+i)->molecule.nDihedrals;j++)
		{
			fgets(buffer, 200, fp);
			sscanf(buffer, "%d %d %d %d", 
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom1),
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom2),
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom3),
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom4),
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->type));
			(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom1)--;
			(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom2)--;
			(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom3)--;
			(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom4)--;
		}

		//Initiate Impropers
		for(j=0;j<(system->moleculeTypes+i)->molecule.nImpropers;j++)
		{
			fgets(buffer, 200, fp);
			sscanf(buffer, "%d %d %d %d", 
								&(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom1),
								&(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom2),
								&(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom3),
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom4),
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->type));
			(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom1)--;
			(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom2)--;
			(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom3)--;
			(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom4)--;
		}
	}
	return 1;
}

int CreateMolecules(SYSTEM *system)
{
	int i=0;
	int j=0;
	int k=0;
	int countPerDim;
	double molecularPosx,molecularPosy,molecularPosz;
	
	system->nAllMolecules = 0;
	for(i=0; i<system->nMoleculeTypes; i++)
	{
		system->nAllMolecules += (system->moleculeTypes+i)->nMolecules;
	}
	
	//Allocate memory for allMolecules
	system->allMolecules = malloc(system->nAllMolecules*sizeof(MOLECULE));
	
	//Initialize allMolecules
	int iAll=0;
	for(i=0; i<system->nMoleculeTypes; i++)
	{
		for(j=0; j<(system->moleculeTypes+i)->nMolecules; j++)
		{
			(system->allMolecules+iAll)->nAtoms = 
					(system->moleculeTypes+i)->molecule.nAtoms;
			(system->allMolecules+iAll)->nBonds = 
					(system->moleculeTypes+i)->molecule.nBonds;
			(system->allMolecules+iAll)->nAngles = 
					(system->moleculeTypes+i)->molecule.nAngles; 
			(system->allMolecules+iAll)->nDihedrals = 
					(system->moleculeTypes+i)->molecule.nDihedrals;
			(system->allMolecules+iAll)->nImpropers = 
					(system->moleculeTypes+i)->molecule.nImpropers;
					
			(system->allMolecules+iAll)->atoms = malloc((system->allMolecules+iAll)->nAtoms*sizeof(ATOM));
			
			//Initiate molecules at random position
			double molecularPosx=(double)rand()/(double)RAND_MAX*system->dimension;
			double molecularPosy=(double)rand()/(double)RAND_MAX*system->dimension;
			double molecularPosz=(double)rand()/(double)RAND_MAX*system->dimension;
			//Initiate molecules at grid position
/*			if(system->is2D)
			{
				countPerDim = ceil(sqrt((system->moleculeTypes+i)->nMolecules));
				molecularPosx = (double)(j%countPerDim)/(double)countPerDim*system->dimension;
				molecularPosy = (double)(j/countPerDim)/(double)countPerDim*system->dimension;
				molecularPosz = 0.0;
			}
			else
			{
				countPerDim = ceil(pow((system->moleculeTypes+i)->nMolecules, 1.0/3.0));
				molecularPosx= (double)(j%(countPerDim*countPerDim)%countPerDim)/(double)countPerDim*system->dimension;
				molecularPosy= (double)(j%(countPerDim*countPerDim)/countPerDim)/(double)countPerDim*system->dimension;
				molecularPosz= (double)(j/(countPerDim*countPerDim))/(double)countPerDim*system->dimension;
			}
*/
			for(k=0; k<(system->allMolecules+iAll)->nAtoms; k++)
			{
				((system->allMolecules+iAll)->atoms+k)->x = 
						((system->moleculeTypes+i)->molecule.atoms+k)->x + molecularPosx; 
				((system->allMolecules+iAll)->atoms+k)->y = 
						((system->moleculeTypes+i)->molecule.atoms+k)->y + molecularPosy; 
				((system->allMolecules+iAll)->atoms+k)->z = 
						((system->moleculeTypes+i)->molecule.atoms+k)->z + molecularPosz;
				((system->allMolecules+iAll)->atoms+k)->oldx = 
						((system->allMolecules+iAll)->atoms+k)->x; 
				((system->allMolecules+iAll)->atoms+k)->oldy = 
						((system->allMolecules+iAll)->atoms+k)->y; 
				((system->allMolecules+iAll)->atoms+k)->oldz = 
						((system->allMolecules+iAll)->atoms+k)->z; 
/*			
				//Init zero velocity	
				((system->allMolecules+iAll)->atoms+k)->vx = 
						((system->moleculeTypes+i)->molecule.atoms+k)->vx ; 
				((system->allMolecules+iAll)->atoms+k)->vy = 
						((system->moleculeTypes+i)->molecule.atoms+k)->vy ; 
				((system->allMolecules+iAll)->atoms+k)->vz = 
						((system->moleculeTypes+i)->molecule.atoms+k)->vz ; 
*/			
				//Init random velocity			
				((system->allMolecules+iAll)->atoms+k)->vx = 
						(double)rand()/(double)RAND_MAX ; 
				((system->allMolecules+iAll)->atoms+k)->vy = 
						(double)rand()/(double)RAND_MAX ; 
				((system->allMolecules+iAll)->atoms+k)->vz = 
						(double)rand()/(double)RAND_MAX ; 
						
				((system->allMolecules+iAll)->atoms+k)->mol = iAll;
				((system->allMolecules+iAll)->atoms+k)->ax = 
						((system->moleculeTypes+i)->molecule.atoms+k)->ax ; 
				((system->allMolecules+iAll)->atoms+k)->ay = 
						((system->moleculeTypes+i)->molecule.atoms+k)->ay ; 
				((system->allMolecules+iAll)->atoms+k)->az = 
						((system->moleculeTypes+i)->molecule.atoms+k)->az ; 
				((system->allMolecules+iAll)->atoms+k)->mass = 
						((system->moleculeTypes+i)->molecule.atoms+k)->mass ; 
				((system->allMolecules+iAll)->atoms+k)->type = 
						((system->moleculeTypes+i)->molecule.atoms+k)->type ; 
				((system->allMolecules+iAll)->atoms+k)->charge = 
						((system->moleculeTypes+i)->molecule.atoms+k)->charge ; 
			}
		
			//Initialize every bonds
			(system->allMolecules+iAll)->bonds = malloc((system->allMolecules+iAll)->nBonds*sizeof(BOND));
			for(k=0; k<(system->allMolecules+iAll)->nBonds; k++)
			{
				((system->allMolecules+iAll)->bonds+k)->type = 
						((system->moleculeTypes+i)->molecule.bonds+k)->type ;
				((system->allMolecules+iAll)->bonds+k)->idxAtom1 = 
						((system->moleculeTypes+i)->molecule.bonds+k)->idxAtom1 ;
				((system->allMolecules+iAll)->bonds+k)->idxAtom2 = 
						((system->moleculeTypes+i)->molecule.bonds+k)->idxAtom2 ;
				((system->allMolecules+iAll)->bonds+k)->atom1 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.bonds+k)->idxAtom1) ;
				((system->allMolecules+iAll)->bonds+k)->atom2 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.bonds+k)->idxAtom2) ;
			}
		
			(system->allMolecules+iAll)->angles = malloc((system->allMolecules+iAll)->nAngles*sizeof(ANGLE));
			for(k=0; k<(system->allMolecules+iAll)->nAngles; k++)
			{
				((system->allMolecules+iAll)->angles+k)->type = 
						((system->moleculeTypes+i)->molecule.angles+k)->type ;
				((system->allMolecules+iAll)->angles+k)->idxAtom1 = 
						((system->moleculeTypes+i)->molecule.angles+k)->idxAtom1 ;
				((system->allMolecules+iAll)->angles+k)->idxAtom2 = 
						((system->moleculeTypes+i)->molecule.angles+k)->idxAtom2 ;
				((system->allMolecules+iAll)->angles+k)->idxAtom3 = 
						((system->moleculeTypes+i)->molecule.angles+k)->idxAtom3 ;
				((system->allMolecules+iAll)->angles+k)->atom1 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.angles+k)->idxAtom1) ;
				((system->allMolecules+iAll)->angles+k)->atom2 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.angles+k)->idxAtom2) ;
				((system->allMolecules+iAll)->angles+k)->atom3 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.angles+k)->idxAtom3) ;
			}

			(system->allMolecules+iAll)->dihedrals = malloc((system->allMolecules+iAll)->nDihedrals*sizeof(DIHEDRAL));
			for(k=0; k<(system->allMolecules+iAll)->nDihedrals; k++)
			{
				((system->allMolecules+iAll)->dihedrals+k)->type = 
						((system->moleculeTypes+i)->molecule.dihedrals+k)->type ;
				((system->allMolecules+iAll)->dihedrals+k)->idxAtom1 = 
						((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom1 ;
				((system->allMolecules+iAll)->dihedrals+k)->idxAtom2 = 
						((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom2 ;
				((system->allMolecules+iAll)->dihedrals+k)->idxAtom3 = 
						((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom3 ;
				((system->allMolecules+iAll)->dihedrals+k)->idxAtom3 = 
						((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom3 ;
				((system->allMolecules+iAll)->dihedrals+k)->atom1 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom1) ;
				((system->allMolecules+iAll)->dihedrals+k)->atom2 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom2) ;
				((system->allMolecules+iAll)->dihedrals+k)->atom3 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom3) ;
				((system->allMolecules+iAll)->dihedrals+k)->atom3 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.dihedrals+k)->idxAtom3) ;
			}

			(system->allMolecules+iAll)->impropers = malloc((system->allMolecules+iAll)->nImpropers*sizeof(IMPROPER));
			for(k=0; k<(system->allMolecules+iAll)->nImpropers; k++)
			{
				((system->allMolecules+iAll)->impropers+k)->type = 
						((system->moleculeTypes+i)->molecule.impropers+k)->type ;
				((system->allMolecules+iAll)->impropers+k)->idxAtom1 = 
						((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom1 ;
				((system->allMolecules+iAll)->impropers+k)->idxAtom2 = 
						((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom2 ;
				((system->allMolecules+iAll)->impropers+k)->idxAtom3 = 
						((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom3 ;
				((system->allMolecules+iAll)->impropers+k)->idxAtom3 = 
						((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom3 ;
				((system->allMolecules+iAll)->impropers+k)->atom1 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom1) ;
				((system->allMolecules+iAll)->impropers+k)->atom2 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom2) ;
				((system->allMolecules+iAll)->impropers+k)->atom3 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom3) ;
				((system->allMolecules+iAll)->impropers+k)->atom3 = 
						(	(system->allMolecules+iAll)->atoms + 
							((system->moleculeTypes+i)->molecule.impropers+k)->idxAtom3) ;
			}
			iAll++;
		}
	}
	return 1;
}

int	CalculateBondForce(SYSTEM *system)
{
	int i=0;
	int j=0;
	double dx,dy,dz,r,dr,forceCoef, K, r0;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nBonds; j++)
		{
			K  =  (system->bondTypes+((system->allMolecules+i)->bonds+j)->type)->K;
			r0 =  (system->bondTypes+((system->allMolecules+i)->bonds+j)->type)->r0;
			dx =  ((system->allMolecules+i)->bonds+j)->atom1->x - 
						((system->allMolecules+i)->bonds+j)->atom2->x;
			dy =  ((system->allMolecules+i)->bonds+j)->atom1->y - 
						((system->allMolecules+i)->bonds+j)->atom2->y;
			dz =  ((system->allMolecules+i)->bonds+j)->atom1->z - 
						((system->allMolecules+i)->bonds+j)->atom2->z;

			//get minimun image
			if (fabs(dx) > system->dimension/2.0) 
			{
				if (dx < 0.0) dx += system->dimension;
				else dx -= system->dimension;
			}
			if (fabs(dy) > system->dimension/2.0) 
			{
				if (dy < 0.0) dy += system->dimension;
				else dy -= system->dimension;
			}
			if (fabs(dz) > system->dimension/2.0) 
			{
				if (dz < 0.0) dz += system->dimension;
				else dz -= system->dimension;
			}

			r = sqrt(dx*dx+dy*dy+dz*dz);
			dr = r -r0;
			forceCoef = 2.0 * K * dr / r *4.184e-4;
			
			system->potentialEnergy += K * dr *dr *4.184e-4;
			
			((system->allMolecules+i)->bonds+j)->atom1->ax -= forceCoef*dx/((system->allMolecules+i)->bonds+j)->atom1->mass;
			((system->allMolecules+i)->bonds+j)->atom1->ay -= forceCoef*dy/((system->allMolecules+i)->bonds+j)->atom1->mass;
			((system->allMolecules+i)->bonds+j)->atom1->az -= forceCoef*dz/((system->allMolecules+i)->bonds+j)->atom1->mass;
			((system->allMolecules+i)->bonds+j)->atom2->ax += forceCoef*dx/((system->allMolecules+i)->bonds+j)->atom2->mass;
			((system->allMolecules+i)->bonds+j)->atom2->ay += forceCoef*dy/((system->allMolecules+i)->bonds+j)->atom2->mass;
			((system->allMolecules+i)->bonds+j)->atom2->az += forceCoef*dz/((system->allMolecules+i)->bonds+j)->atom2->mass;//
		}
	}
	return 1;
}

int CalculateAngleForce(SYSTEM *system)
{
	int i,j;
	double K, theta0, dx1, dy1, dz1, rsq1, r1, dx2 ,dy2, dz2, rsq2, r2, cosTheta, sinTheta, dtheta, tk, a, a11, a12, a22, f1x, f1y, f1z, f3x, f3y, f3z;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAngles; j++)
		{
			K  =  (system->angleTypes+((system->allMolecules+i)->angles+j)->type)->K;
			theta0 =  (system->angleTypes+((system->allMolecules+i)->angles+j)->type)->theta0/180.0*3.14159265357989;
			dx1 =  ((system->allMolecules+i)->angles+j)->atom1->x - 
						((system->allMolecules+i)->angles+j)->atom2->x;
			dy1 =  ((system->allMolecules+i)->angles+j)->atom1->y - 
						((system->allMolecules+i)->angles+j)->atom2->y;
			dz1 =  ((system->allMolecules+i)->angles+j)->atom1->z - 
						((system->allMolecules+i)->angles+j)->atom2->z;

			//get minimun image
			if (fabs(dx1) > system->dimension/2.0) 
			{
				if (dx1 < 0.0) dx1 += system->dimension;
				else dx1 -= system->dimension;
			}
			if (fabs(dy1) > system->dimension/2.0) 
			{
				if (dy1 < 0.0) dy1 += system->dimension;
				else dy1 -= system->dimension;
			}
			if (fabs(dz1) > system->dimension/2.0) 
			{
				if (dz1 < 0.0) dz1 += system->dimension;
				else dz1 -= system->dimension;
			}

			rsq1 = dx1*dx1+dy1*dy1+dz1*dz1;
			r1 = sqrt(rsq1);


			dx2 =  ((system->allMolecules+i)->angles+j)->atom3->x - 
						((system->allMolecules+i)->angles+j)->atom2->x;
			dy2 =  ((system->allMolecules+i)->angles+j)->atom3->y - 
						((system->allMolecules+i)->angles+j)->atom2->y;
			dz2 =  ((system->allMolecules+i)->angles+j)->atom3->z - 
						((system->allMolecules+i)->angles+j)->atom2->z;

			//get minimun image
			if (fabs(dx2) > system->dimension/2.0) 
			{
				if (dx2 < 0.0) dx2 += system->dimension;
				else dx2 -= system->dimension;
			}
			if (fabs(dy2) > system->dimension/2.0) 
			{
				if (dy2 < 0.0) dy2 += system->dimension;
				else dy2 -= system->dimension;
			}
			if (fabs(dz2) > system->dimension/2.0) 
			{
				if (dz2 < 0.0) dz2 += system->dimension;
				else dz2 -= system->dimension;
			}

			rsq2 = dx2*dx2+dy2*dy2+dz2*dz2;
			r2 = sqrt(rsq2);


			cosTheta = dx1*dx2 + dy1*dy2 + dz1*dz2;
			cosTheta /= r1*r2;
			if (cosTheta > 1.0) cosTheta= 1.0;
			if (cosTheta <-1.0) cosTheta=-1.0;
			
			sinTheta = sqrt(1.0 -cosTheta*cosTheta);
			if(sinTheta < 0.001) sinTheta=0.001;
			dtheta = acos(cosTheta) - theta0;
			
			tk = K*dtheta;
			
			a   = -2.0*tk/sinTheta*4.184e-4;
			a11 = a*cosTheta / rsq1;
			a12 = -a /(r1*r2);
			a22 = a*cosTheta / rsq2;
			
			f1x = a11*dx1 + a12*dx2;
			f1y = a11*dy1 + a12*dy2;
			f1z = a11*dz1 + a12*dz2;
			f3x = a22*dx2 + a12*dx1;
			f3y = a22*dy2 + a12*dy1;
			f3z = a22*dz2 + a12*dz1;
			
			system->potentialEnergy += K * dtheta *dtheta *4.184e-4;
			((system->allMolecules+i)->angles+j)->atom1->ax += f1x/((system->allMolecules+i)->angles+j)->atom1->mass;
			((system->allMolecules+i)->angles+j)->atom1->ay += f1y/((system->allMolecules+i)->angles+j)->atom1->mass;
			((system->allMolecules+i)->angles+j)->atom1->az += f1z/((system->allMolecules+i)->angles+j)->atom1->mass;
			((system->allMolecules+i)->angles+j)->atom3->ax += f3x/((system->allMolecules+i)->angles+j)->atom3->mass;
			((system->allMolecules+i)->angles+j)->atom3->ay += f3y/((system->allMolecules+i)->angles+j)->atom3->mass;
			((system->allMolecules+i)->angles+j)->atom3->az += f3z/((system->allMolecules+i)->angles+j)->atom3->mass;
			((system->allMolecules+i)->angles+j)->atom2->ax -= (f1x+f3x)/((system->allMolecules+i)->angles+j)->atom2->mass;
			((system->allMolecules+i)->angles+j)->atom2->ay -= (f1y+f3y)/((system->allMolecules+i)->angles+j)->atom2->mass;
			((system->allMolecules+i)->angles+j)->atom2->az -= (f1z+f3z)/((system->allMolecules+i)->angles+j)->atom2->mass;
		}
	}
	
	return 1;
}

//Pair force implemented grid sckeme
int CalculatePairForce(SYSTEM *system)
{
	int i,j,k,l,m,n,o,p,im,in,io;
	double dx,dy,dz,r2,r,r6,LJCoef,coulCoef,forceCoef;
	double Cepsilon = 0.01;
	double epsilon = 0.001;
	double soft = 0.001;
	ATOM *atom1;
	ATOM *atom2;
	for(i=0; i<system->nSlice; i++)
	{
		for(j=0; j<system->nSlice; j++)
		{
			for(k=0; k<system->nSlice; k++)
			{
				for(l=0; l<system->gridCount[i*system->nSlice*system->nSlice+j*system->nSlice+k]; l++)
				{
					atom1 = system->grid[i*system->nSlice*system->nSlice+j*system->nSlice+k][l];
					for(im=i-1; im<i+2; im++)
					{
						for(in=j-1; in<j+2; in++)
						{
							for(io=k-1; io<k+2; io++)
							{
								m=im;
								n=in;
								o=io;
								if(m<0) m+=system->nSlice;
								if(n<0) n+=system->nSlice;
								if(o<0) o+=system->nSlice;
								if(m>=system->nSlice) m-=system->nSlice;
								if(n>=system->nSlice) n-=system->nSlice;
								if(o>=system->nSlice) o-=system->nSlice;
								for(p=0; p<system->gridCount[m*system->nSlice*system->nSlice+n*system->nSlice+o]; p++)
								{
									atom2 = system->grid[m*system->nSlice*system->nSlice+n*system->nSlice+o][p];
									if(atom1->mol == atom2->mol) continue;
									dx = atom2->x - atom1->x;
									dy = atom2->y - atom1->y;
									dz = atom2->z - atom1->z;
									
									//get minimun image
									if (fabs(dx) > system->dimension/2.0) 
									{
										if (dx < 0.0) dx += system->dimension;
										else dx -= system->dimension;
									}
									if (fabs(dy) > system->dimension/2.0) 
									{
										if (dy < 0.0) dy += system->dimension;
										else dy -= system->dimension;
									}
									if (fabs(dz) > system->dimension/2.0) 
									{
										if (dz < 0.0) dz += system->dimension;
										else dz -= system->dimension;
									}

									r2 = dx*dx + dy*dy + dz*dz + soft;
									if (r2 > system->rCut*system->rCut) continue;
									r = sqrt(r2);

									coulCoef = Cepsilon*atom2->charge*atom1->charge/r/r2;

									r6 = r2*r2*r2;
									
									LJCoef = 4.0 * epsilon * (6.0 - 12.0/r6)/r6/r2;
						
									forceCoef = (LJCoef-coulCoef) / atom1->mass;
//									if (forceCoef < -0.01) forceCoef=-0.01;
									atom1->ax += forceCoef * dx;
									atom1->ay += forceCoef * dy;
									atom1->az += forceCoef * dz;
									
									system->potentialEnergy -= 4.0 * epsilon * (1.0 - 1.0/r6)/r6;
									system->potentialEnergy += Cepsilon*atom2->charge*atom1->charge/r;
									
								}
							}
						}
					}
				}
			}
		}
	}
	return 1;
}

//Global pair force 
int	CalculatePairForceOld(SYSTEM *system)
{
	int i=0;
	int j=0;
	int k=0;
	int l=0;
	double dx,dy,dz,r2,r4,r8,r14,forceCoef;
	double epsilon = 0.001;
	double soft = 0.001;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			for(k=0; k<system->nAllMolecules; k++)
			{
				for(l=0; l<(system->allMolecules+i)->nAtoms; l++)
				{
					if(i == k) continue;
					dx = 	((system->allMolecules+k)->atoms+l)->x - 
								((system->allMolecules+i)->atoms+j)->x;
					dy = 	((system->allMolecules+k)->atoms+l)->y - 
								((system->allMolecules+i)->atoms+j)->y;
					dz = 	((system->allMolecules+k)->atoms+l)->z - 
								((system->allMolecules+i)->atoms+j)->z;

					//get minimun image
					if (fabs(dx) > system->dimension/2.0) 
					{
						if (dx < 0.0) dx += system->dimension;
						else dx -= system->dimension;
					}
					if (fabs(dy) > system->dimension/2.0) 
					{
						if (dy < 0.0) dy += system->dimension;
						else dy -= system->dimension;
					}
					if (fabs(dz) > system->dimension/2.0) 
					{
						if (dz < 0.0) dz += system->dimension;
						else dz -= system->dimension;
					}

					r2 = dx*dx + dy*dy + dz*dz + soft;
					r4 = r2*r2;
					r8 = r4*r4;
					r14 = r8*r4*r2;
					forceCoef = 4.0 * epsilon * (1.0/r8 - 1.0/r14) / ((system->allMolecules+i)->atoms+j)->mass;
					((system->allMolecules+i)->atoms+j)->ax -= forceCoef * dx;
					((system->allMolecules+i)->atoms+j)->ay -= forceCoef * dy;
					((system->allMolecules+i)->atoms+j)->az -= forceCoef * dz;
					system->potentialEnergy -= 4.0 * epsilon * (1.0 - 1.0/r4/r2)/r4/r2;
					
				}
			}
		}
	}
	return 1;
}
int ReleaseGridList(SYSTEM *system)
{
	int i;
	for(i=0; i<system->nGrids; i++)
	{
		free(system->grid[i]);
	}
	free(system->grid);
	free(system->gridCount);
	return 1;
}
int CreateGridList(SYSTEM *system)
{
	int i,j,k,ix,iy,iz,idx;
	
	system->nSlice = floor(system->dimension / system->rCut);
	system->nGrids = system->nSlice*system->nSlice*system->nSlice;
	system->grid = malloc(system->nGrids*sizeof(ATOM**));
	system->gridCount = malloc(system->nGrids*sizeof(int));
	
	//get grid counts
	memset(system->gridCount, 0, system->nGrids*sizeof(int));
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			ix = floor(((system->allMolecules+i)->atoms+j)->x /
				(system->dimension/(double)(system->nSlice)));
			iy = floor(((system->allMolecules+i)->atoms+j)->y /
				(system->dimension/(double)(system->nSlice)));
			iz = floor(((system->allMolecules+i)->atoms+j)->z /
				(system->dimension/(double)(system->nSlice)));
			if(ix<0) ix=0;
			if(iy<0) iy=0;
			if(iz<0) iz=0;
			if(ix>=system->nSlice) ix=system->nSlice-1;
			if(iy>=system->nSlice) iy=system->nSlice-1;
			if(iz>=system->nSlice) iz=system->nSlice-1;
			
			system->gridCount[iz*system->nSlice*system->nSlice+iy*system->nSlice+ix]++;
		}
	}
	
	//allocate grid space
	for(i=0; i<system->nSlice; i++)
	{
		for(j=0; j<system->nSlice; j++)
		{
			for(k=0; k<system->nSlice; k++)
			{
				idx = i*system->nSlice*system->nSlice+j*system->nSlice+k;
				system->grid[idx] = 
					malloc(system->gridCount[idx]*sizeof(ATOM**));
			}
		}
	}

	//set grid
	memset(system->gridCount, 0, system->nGrids*sizeof(int));
	
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			ix = floor(((system->allMolecules+i)->atoms+j)->x /
				(system->dimension/(double)(system->nSlice)));
			iy = floor(((system->allMolecules+i)->atoms+j)->y /
				(system->dimension/(double)(system->nSlice)));
			iz = floor(((system->allMolecules+i)->atoms+j)->z /
				(system->dimension/(double)(system->nSlice)));
			if(ix<0) ix=0;
			if(iy<0) iy=0;
			if(iz<0) iz=0;
			if(ix>=system->nSlice) ix=system->nSlice-1;
			if(iy>=system->nSlice) iy=system->nSlice-1;
			if(iz>=system->nSlice) iz=system->nSlice-1;
			idx = iz*system->nSlice*system->nSlice+iy*system->nSlice+ix;
//			printf("%d %d %d %d %d %lf\n",ix,iy,iz,idx, system->gridCount[idx], ((system->allMolecules+i)->atoms+j)->z);
			system->grid[idx][system->gridCount[idx]]=(system->allMolecules+i)->atoms+j;
			system->gridCount[idx]++;
		}
	}
/*	
	for(i=0; i<system->nSlice; i++)
	{
		for(j=0; j<system->nSlice; j++)
		{
			for(k=0; k<system->nSlice; k++)
			{
				printf("%d ",system->gridCount[i*system->nSlice*system->nSlice+j*system->nSlice+k]);
			}
			printf("\n");
		}
		printf("\n");
	}
*/	
	return 1;
}
int	CalculateForce(SYSTEM *system)
{
	int i=0;
	int j=0;
	//Reset Accelerations
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			((system->allMolecules+i)->atoms+j)->ax = 0.0;
			((system->allMolecules+i)->atoms+j)->ay = 0.0;
			((system->allMolecules+i)->atoms+j)->az = 0.0;
		}
	}
	system->potentialEnergy = 0.0;
	CalculateBondForce(system);
	CalculateAngleForce(system);
	CalculatePairForce(system);
//	getchar();
	return 1;
}


int	Integrate(SYSTEM *system)
{
	int i=0;
	int j=0;
	double newx,newy,newz,dx,dy,dz;
	double dt=0.5;
	ATOM *atom;
//	double max = 0.0;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			atom = (system->allMolecules+i)->atoms+j;
			newx = 2.0 *	atom->x - atom->oldx + atom->ax * dt * dt;
			newy = 2.0 *	atom->y - atom->oldy + atom->ay * dt * dt;
			newz = 2.0 *	atom->z - atom->oldz + atom->az * dt * dt;
			
			atom->oldx = atom->x;
			atom->oldy = atom->y;
			atom->oldz = atom->z;
			dx = newx - atom->oldx;
			dy = newy - atom->oldy;
			dz = newz - atom->oldz;
			//get minimun image
			if (fabs(dx) > system->dimension/2.0) 
			{
				if (dx < 0.0) dx += system->dimension;
				else dx -= system->dimension;
			}
			if (fabs(dy) > system->dimension/2.0) 
			{
				if (dy < 0.0) dy += system->dimension;
				else dy -= system->dimension;
			}
			if (fabs(dz) > system->dimension/2.0) 
			{
				if (dz < 0.0) dz += system->dimension;
				else dz -= system->dimension;
			}

			atom->vx = dx/dt;
			atom->vy = dy/dt;
			atom->vz = dz/dt;
//			if(fabs(atom->ax) > 0.005) atom->ax=0.0;
//			if(fabs(atom->ay) > 0.005) atom->ay=0.0;
//			if(fabs(atom->az) > 0.005) atom->az=0.0;
										
			
			if(newx<0) newx+=system->dimension;
			if(newy<0) newy+=system->dimension;
			if(newz<0) newz+=system->dimension;
			if(newx>=system->dimension) newx-=system->dimension;
			if(newy>=system->dimension) newy-=system->dimension;
			if(newz>=system->dimension) newz-=system->dimension;
			atom->x = newx;
			atom->y = newy;
			atom->z = newz;
		}
	}
	return 1;
}

int	CalculateKineticEnergy(SYSTEM *system)
{
	int i=0;
	int j=0;
	
	system->kineticEnergy = 0.0;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			system->kineticEnergy += 	((system->allMolecules+i)->atoms+j)->vx*
																((system->allMolecules+i)->atoms+j)->vx*
																((system->allMolecules+i)->atoms+j)->mass*0.5;
			system->kineticEnergy += 	((system->allMolecules+i)->atoms+j)->vy*
																((system->allMolecules+i)->atoms+j)->vy*
																((system->allMolecules+i)->atoms+j)->mass*0.5;
			system->kineticEnergy += 	((system->allMolecules+i)->atoms+j)->vz*
																((system->allMolecules+i)->atoms+j)->vz*
																((system->allMolecules+i)->atoms+j)->mass*0.5;
		}
	}
	return 1;
}

int	UpdateMolecules(SYSTEM *system)
{
	Integrate(system);
	CreateGridList(system);
	CalculateForce(system);
//	RenderMolecules(system);
	ReleaseGridList(system);
	CalculateKineticEnergy(system);
	return 1;
}

int	IntegrateRelaxation(SYSTEM *system)
{
	int i=0;
	int j=0;
	double newx,newy,newz,dx,dy,dz;
	double dt=0.5;
	ATOM *atom;
//	double max = 0.0;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			atom = (system->allMolecules+i)->atoms+j;
			
			atom->ax -= atom->vx;
			atom->ay -= atom->vy;
			atom->az -= atom->vz;
			
			newx = 2.0 *	atom->x - atom->oldx + atom->ax * dt * dt;
			newy = 2.0 *	atom->y - atom->oldy + atom->ay * dt * dt;
			newz = 2.0 *	atom->z - atom->oldz + atom->az * dt * dt;
			
			atom->oldx = atom->x;
			atom->oldy = atom->y;
			atom->oldz = atom->z;
			dx = newx - atom->oldx*0.999;
			dy = newy - atom->oldy*0.999;
			dz = newz - atom->oldz*0.999;
			//get minimun image
			if (fabs(dx) > system->dimension/2.0) 
			{
				if (dx < 0.0) dx += system->dimension;
				else dx -= system->dimension;
			}
			if (fabs(dy) > system->dimension/2.0) 
			{
				if (dy < 0.0) dy += system->dimension;
				else dy -= system->dimension;
			}
			if (fabs(dz) > system->dimension/2.0) 
			{
				if (dz < 0.0) dz += system->dimension;
				else dz -= system->dimension;
			}

			atom->vx = dx/dt;
			atom->vy = dy/dt;
			atom->vz = dz/dt;
//			if(fabs(atom->ax) > 0.005) atom->ax=0.0;
//			if(fabs(atom->ay) > 0.005) atom->ay=0.0;
//			if(fabs(atom->az) > 0.005) atom->az=0.0;
										
			
			if(newx<0) newx+=system->dimension;
			if(newy<0) newy+=system->dimension;
			if(newz<0) newz+=system->dimension;
			if(newx>=system->dimension) newx-=system->dimension;
			if(newy>=system->dimension) newy-=system->dimension;
			if(newz>=system->dimension) newz-=system->dimension;
			atom->x = newx;
			atom->y = newy;
			atom->z = newz;
		}
	}
	return 1;
}


int RelaxMolecules(SYSTEM *system)
{
	IntegrateRelaxation(system);
	CreateGridList(system);
	CalculateForce(system);
//	RenderMolecules(system);
	ReleaseGridList(system);
	CalculateKineticEnergy(system);
	return 1;
}

int ResetVelocity(SYSTEM *system)
{
	int i,j;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			((system->allMolecules+i)->atoms+j)->vx=0.0;
			((system->allMolecules+i)->atoms+j)->vy=0.0;
			((system->allMolecules+i)->atoms+j)->vz=0.0;
			((system->allMolecules+i)->atoms+j)->ax=0.0;
			((system->allMolecules+i)->atoms+j)->ay=0.0;
			((system->allMolecules+i)->atoms+j)->az=0.0;
			((system->allMolecules+i)->atoms+j)->oldx=((system->allMolecules+i)->atoms+j)->x;
			((system->allMolecules+i)->atoms+j)->oldy=((system->allMolecules+i)->atoms+j)->y;
			((system->allMolecules+i)->atoms+j)->oldz=((system->allMolecules+i)->atoms+j)->z;
		}
	}
	return 1;
}

int ReleaseMolecules(SYSTEM *system) //FIXME
{
	int i=0;
	
	//Release contents of every molecule
	for(i=0; i<system->nAllMolecules; i++)
	{
		free((system->allMolecules+i)->atoms);
		free((system->allMolecules+i)->bonds);
		free((system->allMolecules+i)->angles);
		free((system->allMolecules+i)->dihedrals);
		free((system->allMolecules+i)->impropers);
		(system->allMolecules+i)->atoms=NULL;
		(system->allMolecules+i)->bonds=NULL;
		(system->allMolecules+i)->angles=NULL;
		(system->allMolecules+i)->dihedrals=NULL;
		(system->allMolecules+i)->impropers=NULL;
	}
	
	//Release allMolecules
	free(system->allMolecules);
	system->allMolecules = NULL;
	return 1;
}

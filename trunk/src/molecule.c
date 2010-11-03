#include "molecule.h"
#define ATOM_TYPE_MAX 21

int InitSystem(SYSTEM *system, char *infile)
{
	int i=0;
	int j=0;
	int k=0;
	
	const char AtomType[ATOM_TYPE_MAX][4]=
	{
		"X",
		"H","He","Li","Be","B","C","N","O","F","Ne",
		"Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
	};
	FILE *fp = fopen(infile,"r");
	char buffer[200];
	char type[4];
	int idummy;
	system->r2min = 9e99; //DEBUG FIXME
	system->dimension = 20.0;
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
			sscanf(buffer, "%lf %lf %lf %s %lf", 
								&(((system->moleculeTypes+i)->molecule.atoms+j)->x),
								&(((system->moleculeTypes+i)->molecule.atoms+j)->y),
								&(((system->moleculeTypes+i)->molecule.atoms+j)->z),
								type,
								&(((system->moleculeTypes+i)->molecule.atoms+j)->mass));
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
			((system->moleculeTypes+i)->molecule.atoms+j)->charge = 0.0;
			
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
			sscanf(buffer, "%d %d %d %d", 
								&(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom1),
								&(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom2),
								&(((system->moleculeTypes+i)->molecule.bonds+j)->type),
								&idummy);
			(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom1)--;
			(((system->moleculeTypes+i)->molecule.bonds+j)->idxAtom2)--;
		}

		//Initiate Angles
		for(j=0;j<(system->moleculeTypes+i)->molecule.nAngles;j++)
		{
			fgets(buffer, 200, fp);
			sscanf(buffer, "%d %d %d", 
								&(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom1),
								&(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom2),
								&(((system->moleculeTypes+i)->molecule.angles+j)->idxAtom3));
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
								&(((system->moleculeTypes+i)->molecule.dihedrals+j)->idxAtom4));
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
								&(((system->moleculeTypes+i)->molecule.impropers+j)->idxAtom4));
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
/*			
			//Initiate molecules at random position
			double molecularPosx=(double)rand()/(double)RAND_MAX*system->dimension;
			double molecularPosy=(double)rand()/(double)RAND_MAX*system->dimension;
			double molecularPosz=(double)rand()/(double)RAND_MAX*system->dimension;
*/			
			//Initiate molecules at grid position
			int countPerDim = ceil(pow((system->moleculeTypes+i)->nMolecules, 1.0/3.0));
			double molecularPosx= (double)(j%(countPerDim*countPerDim)%countPerDim)/(double)countPerDim*system->dimension;
			double molecularPosy= (double)(j%(countPerDim*countPerDim)/countPerDim)/(double)countPerDim*system->dimension;
			double molecularPosz= (double)(j/(countPerDim*countPerDim))/(double)countPerDim*system->dimension;
			
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
	double dx,dy,dz,r,dr,forceCoef;
	double K = 80.0;
	double r0 = 1.15;
	system->potentialEnergy = 0.0;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nBonds; j++)
		{
			dx =  ((system->allMolecules+i)->bonds+j)->atom1->x - 
						((system->allMolecules+i)->bonds+j)->atom2->x;
			dy =  ((system->allMolecules+i)->bonds+j)->atom1->y - 
						((system->allMolecules+i)->bonds+j)->atom2->y;
			dz =  ((system->allMolecules+i)->bonds+j)->atom1->z - 
						((system->allMolecules+i)->bonds+j)->atom2->z;
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

int	CalculatePairForce(SYSTEM *system)
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
					r2 = dx*dx + dy*dy + dz*dz + soft;
					if (r2 < 1.0) r2 = 1.0;  //FIXME
					if(r2 < system->r2min) system->r2min = r2; //DEBUG FIXME
					r4 = r2*r2;
					r8 = r4*r4;
					r14 = r8*r4*r2;
					forceCoef = 4.0 * epsilon * (6.0/r8 - 12.0/r14) / ((system->allMolecules+i)->atoms+j)->mass;
					((system->allMolecules+i)->atoms+j)->ax += forceCoef * dx;
					((system->allMolecules+i)->atoms+j)->ay += forceCoef * dy;
					((system->allMolecules+i)->atoms+j)->az += forceCoef * dz;
					
				}
			}
		}
	}
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
	CalculateBondForce(system);
	CalculatePairForce(system);
	return 1;
}

int	Integrate(SYSTEM *system)
{
	int i=0;
	int j=0;
	double newx,newy,newz;
	double dt=2.0;
	
//	double max = 0.0;
	for(i=0; i<system->nAllMolecules; i++)
	{
		for(j=0; j<(system->allMolecules+i)->nAtoms; j++)
		{
			newx = 2.0 *	((system->allMolecules+i)->atoms+j)->x - 
										((system->allMolecules+i)->atoms+j)->oldx +
										((system->allMolecules+i)->atoms+j)->ax * dt * dt;
			newy = 2.0 *	((system->allMolecules+i)->atoms+j)->y - 
										((system->allMolecules+i)->atoms+j)->oldy +
										((system->allMolecules+i)->atoms+j)->ay * dt * dt;
			newz = 2.0 *	((system->allMolecules+i)->atoms+j)->z - 
										((system->allMolecules+i)->atoms+j)->oldz +
										((system->allMolecules+i)->atoms+j)->az * dt * dt;
/*										
			//Calculate maximum accelaretion
			max = (((system->allMolecules+i)->atoms+j)->ax > max)?((system->allMolecules+i)->atoms+j)->ax:max;
			max = (((system->allMolecules+i)->atoms+j)->ay > max)?((system->allMolecules+i)->atoms+j)->ay:max;
			max = (((system->allMolecules+i)->atoms+j)->az > max)?((system->allMolecules+i)->atoms+j)->az:max;
*/			
			((system->allMolecules+i)->atoms+j)->oldx = ((system->allMolecules+i)->atoms+j)->x;
			((system->allMolecules+i)->atoms+j)->oldy = ((system->allMolecules+i)->atoms+j)->y;
			((system->allMolecules+i)->atoms+j)->oldz = ((system->allMolecules+i)->atoms+j)->z;
			((system->allMolecules+i)->atoms+j)->x = newx;
			((system->allMolecules+i)->atoms+j)->y = newy;
			((system->allMolecules+i)->atoms+j)->z = newz;
			((system->allMolecules+i)->atoms+j)->vx =
				(newx - ((system->allMolecules+i)->atoms+j)->oldx)/dt;
			((system->allMolecules+i)->atoms+j)->vy =
				(newy - ((system->allMolecules+i)->atoms+j)->oldy)/dt;
			((system->allMolecules+i)->atoms+j)->vz =
				(newz - ((system->allMolecules+i)->atoms+j)->oldz)/dt;
		}
	}
//	printf("%lf\n",max);
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
	CalculateForce(system);
	CalculateKineticEnergy(system);
	Integrate(system);
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

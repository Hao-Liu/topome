/*
 * /__ __X  _ \/  __\/  _ \/ \__/|/  __/ Copyright (C) 2010 - Liu Hao
 *   / \ | / \||  \/|| / \|| |\/|||  \  
 *   | | | \_/||  __/| \_/|| |  |||  /_  molecule.c
 *   \_/ \____/\_/   \____/\_/  \|\____\
 * 
 * manipulate molecules
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


#include <libintl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "molecule.h"
#include "system.h"

void get_minimun_image (double *dx, double *dy, double *dz, System *tpm_system)
{
	if (fabs(*dx) > tpm_system->half_dimension) 
	{
		if (*dx < 0.0)
		{
		  *dx += tpm_system->dimension;
		}
		else
		{
		  *dx -= tpm_system->dimension;
		}
	}
	if (fabs(*dy) > tpm_system->half_dimension) 
	{
		if (*dy < 0.0)
		{
		  *dy += tpm_system->dimension;
		}
		else
		{
		  *dy -= tpm_system->dimension;
		}
	}
	if (fabs(*dz) > tpm_system->half_dimension) 
	{
		if (*dz < 0.0)
		{
		  *dz += tpm_system->dimension;
		}
		else
		{
		  *dz -= tpm_system->dimension;
		}
	}
}


int	CalculateBondForce(System *tpm_system)
{
	int i=0;
	int j=0;
	double dx,dy,dz,r,dr,forceCoef, K, r0;
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nBonds; j++)
		{
			K  =  (tpm_system->bond_type+((tpm_system->molecule+i)->bonds+j)->type)->K;
			r0 =  (tpm_system->bond_type+((tpm_system->molecule+i)->bonds+j)->type)->r0;
			dx =  ((tpm_system->molecule+i)->bonds+j)->atom1->x - 
						((tpm_system->molecule+i)->bonds+j)->atom2->x;
			dy =  ((tpm_system->molecule+i)->bonds+j)->atom1->y - 
						((tpm_system->molecule+i)->bonds+j)->atom2->y;
			dz =  ((tpm_system->molecule+i)->bonds+j)->atom1->z - 
						((tpm_system->molecule+i)->bonds+j)->atom2->z;
      
      get_minimun_image (&dx, &dy, &dz, tpm_system);

			r = sqrt(dx*dx+dy*dy+dz*dz);
			dr = r -r0;
			forceCoef = 2.0 * K * dr / r *4.184e-4;
			
			tpm_system->potential_energy += K * dr *dr *4.184e-4;
			
			((tpm_system->molecule+i)->bonds+j)->atom1->ax -= forceCoef*dx/((tpm_system->molecule+i)->bonds+j)->atom1->mass;
			((tpm_system->molecule+i)->bonds+j)->atom1->ay -= forceCoef*dy/((tpm_system->molecule+i)->bonds+j)->atom1->mass;
			((tpm_system->molecule+i)->bonds+j)->atom1->az -= forceCoef*dz/((tpm_system->molecule+i)->bonds+j)->atom1->mass;
			((tpm_system->molecule+i)->bonds+j)->atom2->ax += forceCoef*dx/((tpm_system->molecule+i)->bonds+j)->atom2->mass;
			((tpm_system->molecule+i)->bonds+j)->atom2->ay += forceCoef*dy/((tpm_system->molecule+i)->bonds+j)->atom2->mass;
			((tpm_system->molecule+i)->bonds+j)->atom2->az += forceCoef*dz/((tpm_system->molecule+i)->bonds+j)->atom2->mass;//
		}
	}
	return 1;
}

int CalculateAngleForce(System *tpm_system)
{
	int i,j;
	double K, theta0, dx1, dy1, dz1, rsq1, r1, dx2 ,dy2, dz2, rsq2, r2, cosTheta, sinTheta, dtheta, tk, a, a11, a12, a22, f1x, f1y, f1z, f3x, f3y, f3z;
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAngles; j++)
		{
			K  =  (tpm_system->angle_type+((tpm_system->molecule+i)->angles+j)->type)->K;
			theta0 =  (tpm_system->angle_type+((tpm_system->molecule+i)->angles+j)->type)->theta0/180.0*3.14159265357989;
			dx1 =  ((tpm_system->molecule+i)->angles+j)->atom1->x - 
						((tpm_system->molecule+i)->angles+j)->atom2->x;
			dy1 =  ((tpm_system->molecule+i)->angles+j)->atom1->y - 
						((tpm_system->molecule+i)->angles+j)->atom2->y;
			dz1 =  ((tpm_system->molecule+i)->angles+j)->atom1->z - 
						((tpm_system->molecule+i)->angles+j)->atom2->z;

      get_minimun_image (&dx1, &dy1, &dz1, tpm_system);

			rsq1 = dx1*dx1+dy1*dy1+dz1*dz1;
			r1 = sqrt(rsq1);


			dx2 =  ((tpm_system->molecule+i)->angles+j)->atom3->x - 
						((tpm_system->molecule+i)->angles+j)->atom2->x;
			dy2 =  ((tpm_system->molecule+i)->angles+j)->atom3->y - 
						((tpm_system->molecule+i)->angles+j)->atom2->y;
			dz2 =  ((tpm_system->molecule+i)->angles+j)->atom3->z - 
						((tpm_system->molecule+i)->angles+j)->atom2->z;

      get_minimun_image (&dx2, &dy2, &dz2, tpm_system);

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
			
			tpm_system->potential_energy += K * dtheta *dtheta *4.184e-4;
			((tpm_system->molecule+i)->angles+j)->atom1->ax += f1x/((tpm_system->molecule+i)->angles+j)->atom1->mass;
			((tpm_system->molecule+i)->angles+j)->atom1->ay += f1y/((tpm_system->molecule+i)->angles+j)->atom1->mass;
			((tpm_system->molecule+i)->angles+j)->atom1->az += f1z/((tpm_system->molecule+i)->angles+j)->atom1->mass;
			((tpm_system->molecule+i)->angles+j)->atom3->ax += f3x/((tpm_system->molecule+i)->angles+j)->atom3->mass;
			((tpm_system->molecule+i)->angles+j)->atom3->ay += f3y/((tpm_system->molecule+i)->angles+j)->atom3->mass;
			((tpm_system->molecule+i)->angles+j)->atom3->az += f3z/((tpm_system->molecule+i)->angles+j)->atom3->mass;
			((tpm_system->molecule+i)->angles+j)->atom2->ax -= (f1x+f3x)/((tpm_system->molecule+i)->angles+j)->atom2->mass;
			((tpm_system->molecule+i)->angles+j)->atom2->ay -= (f1y+f3y)/((tpm_system->molecule+i)->angles+j)->atom2->mass;
			((tpm_system->molecule+i)->angles+j)->atom2->az -= (f1z+f3z)/((tpm_system->molecule+i)->angles+j)->atom2->mass;
		}
	}
	
	return 1;
}

//Pair force implemented grid scheme
int CalculatePairForce(System *tpm_system)
{
	int i,j,k,l,m,n,o,p,im,in,io;
	double dx,dy,dz,r2,r,r6,LJCoef,coulCoef,forceCoef;
	double Cepsilon = 0.01;
	double epsilon = 0.001;
	double soft = 0.001;
	ATOM *atom1;
	ATOM *atom2;
	for(i=0; i<tpm_system->number_slice; i++)
	{
		for(j=0; j<tpm_system->number_slice; j++)
		{
			for(k=0; k<tpm_system->number_slice; k++)
			{
				for(l=0; l<tpm_system->number_atom_in_grid[i*tpm_system->number_slice*tpm_system->number_slice+j*tpm_system->number_slice+k]; l++)
				{
					atom1 = tpm_system->grid[i*tpm_system->number_slice*tpm_system->number_slice+j*tpm_system->number_slice+k][l];
					for(im=i-1; im<i+2; im++)
					{
						for(in=j-1; in<j+2; in++)
						{
							for(io=k-1; io<k+2; io++)
							{
								m=im;
								n=in;
								o=io;
								if(m<0) m+=tpm_system->number_slice;
								if(n<0) n+=tpm_system->number_slice;
								if(o<0) o+=tpm_system->number_slice;
								if(m>=tpm_system->number_slice) m-=tpm_system->number_slice;
								if(n>=tpm_system->number_slice) n-=tpm_system->number_slice;
								if(o>=tpm_system->number_slice) o-=tpm_system->number_slice;
								for(p=0; p<tpm_system->number_atom_in_grid[m*tpm_system->number_slice*tpm_system->number_slice+n*tpm_system->number_slice+o]; p++)
								{
									atom2 = tpm_system->grid[m*tpm_system->number_slice*tpm_system->number_slice+n*tpm_system->number_slice+o][p];
									if(atom1->mol == atom2->mol) continue;
									dx = atom2->x - atom1->x;
									dy = atom2->y - atom1->y;
									dz = atom2->z - atom1->z;
									
                  get_minimun_image (&dx, &dy, &dz, tpm_system);

									r2 = dx*dx + dy*dy + dz*dz + soft;
									if (r2 > tpm_system->radius_cut*tpm_system->radius_cut) continue;
									r = sqrt(r2);

									coulCoef = Cepsilon*atom2->charge*atom1->charge/r/r2;

									r6 = r2*r2*r2;
									
									LJCoef = 4.0 * epsilon * (6.0 - 12.0/r6)/r6/r2;
						
									forceCoef = (LJCoef-coulCoef) / atom1->mass;
//									if (forceCoef < -0.01) forceCoef=-0.01;
									atom1->ax += forceCoef * dx;
									atom1->ay += forceCoef * dy;
									atom1->az += forceCoef * dz;
									
									tpm_system->potential_energy -= 4.0 * epsilon * (1.0 - 1.0/r6)/r6;
									tpm_system->potential_energy += Cepsilon*atom2->charge*atom1->charge/r;
									
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
int	CalculatePairForceOld(System *tpm_system)
{
	int i=0;
	int j=0;
	int k=0;
	int l=0;
	double dx,dy,dz,r2,r4,r8,r14,forceCoef;
	double epsilon = 0.001;
	double soft = 0.001;
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			for(k=0; k<tpm_system->number_molecule; k++)
			{
				for(l=0; l<(tpm_system->molecule+i)->nAtoms; l++)
				{
					if(i == k) continue;
					dx = 	((tpm_system->molecule+k)->atoms+l)->x - 
								((tpm_system->molecule+i)->atoms+j)->x;
					dy = 	((tpm_system->molecule+k)->atoms+l)->y - 
								((tpm_system->molecule+i)->atoms+j)->y;
					dz = 	((tpm_system->molecule+k)->atoms+l)->z - 
								((tpm_system->molecule+i)->atoms+j)->z;

          get_minimun_image (&dx, &dy, &dz, tpm_system);

					r2 = dx*dx + dy*dy + dz*dz + soft;
					r4 = r2*r2;
					r8 = r4*r4;
					r14 = r8*r4*r2;
					forceCoef = 4.0 * epsilon * (1.0/r8 - 1.0/r14) / ((tpm_system->molecule+i)->atoms+j)->mass;
					((tpm_system->molecule+i)->atoms+j)->ax -= forceCoef * dx;
					((tpm_system->molecule+i)->atoms+j)->ay -= forceCoef * dy;
					((tpm_system->molecule+i)->atoms+j)->az -= forceCoef * dz;
					tpm_system->potential_energy -= 4.0 * epsilon * (1.0 - 1.0/r4/r2)/r4/r2;
					
				}
			}
		}
	}
	return 1;
}

int release_grid_list(System *tpm_system)
{
	int i;
	for(i=0; i<tpm_system->number_grid; i++)
	{
		free(tpm_system->grid[i]);
	}
	free(tpm_system->grid);
	free(tpm_system->number_atom_in_grid);
	return 1;
}

int create_grid_list(System *tpm_system)
{
	int i,j,k,ix,iy,iz,idx;
	
	tpm_system->number_slice = (int)floor(tpm_system->dimension / tpm_system->radius_cut);
	tpm_system->number_grid = tpm_system->number_slice*tpm_system->number_slice*tpm_system->number_slice;
	tpm_system->grid = malloc((long unsigned int)tpm_system->number_grid*sizeof(ATOM**));
	tpm_system->number_atom_in_grid = malloc((long unsigned int)tpm_system->number_grid*sizeof(int));
	
	//get grid counts
	memset(tpm_system->number_atom_in_grid, 0, (long unsigned int)tpm_system->number_grid*sizeof(int));
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			ix = (int)floor(((tpm_system->molecule+i)->atoms+j)->x /
				(tpm_system->dimension/(double)(tpm_system->number_slice)));
			iy = (int)floor(((tpm_system->molecule+i)->atoms+j)->y /
				(tpm_system->dimension/(double)(tpm_system->number_slice)));
			iz = (int)floor(((tpm_system->molecule+i)->atoms+j)->z /
				(tpm_system->dimension/(double)(tpm_system->number_slice)));
			if(ix<0) ix=0;
			if(iy<0) iy=0;
			if(iz<0) iz=0;
			if(ix>=tpm_system->number_slice) ix=tpm_system->number_slice-1;
			if(iy>=tpm_system->number_slice) iy=tpm_system->number_slice-1;
			if(iz>=tpm_system->number_slice) iz=tpm_system->number_slice-1;
			
			tpm_system->number_atom_in_grid[iz*tpm_system->number_slice*tpm_system->number_slice+iy*tpm_system->number_slice+ix]++;
		}
	}
	
	//allocate grid space
	for(i=0; i<tpm_system->number_slice; i++)
	{
		for(j=0; j<tpm_system->number_slice; j++)
		{
			for(k=0; k<tpm_system->number_slice; k++)
			{
				idx = i*tpm_system->number_slice*tpm_system->number_slice+j*tpm_system->number_slice+k;
				tpm_system->grid[idx] = 
					malloc((long unsigned int)tpm_system->number_atom_in_grid[idx]*sizeof(ATOM**));
			}
		}
	}

	//set grid
	memset(tpm_system->number_atom_in_grid, 0, (long unsigned int)tpm_system->number_grid*sizeof(int));
	
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			ix = (int)floor(((tpm_system->molecule+i)->atoms+j)->x /
				(tpm_system->dimension/(double)(tpm_system->number_slice)));
			iy = (int)floor(((tpm_system->molecule+i)->atoms+j)->y /
				(tpm_system->dimension/(double)(tpm_system->number_slice)));
			iz = (int)floor(((tpm_system->molecule+i)->atoms+j)->z /
				(tpm_system->dimension/(double)(tpm_system->number_slice)));
			if(ix<0) ix=0;
			if(iy<0) iy=0;
			if(iz<0) iz=0;
			if(ix>=tpm_system->number_slice) ix=tpm_system->number_slice-1;
			if(iy>=tpm_system->number_slice) iy=tpm_system->number_slice-1;
			if(iz>=tpm_system->number_slice) iz=tpm_system->number_slice-1;
			idx = iz*tpm_system->number_slice*tpm_system->number_slice+iy*tpm_system->number_slice+ix;
//			printf("%d %d %d %d %d %lf\n",ix,iy,iz,idx, tpm_system->number_atom_in_grid[idx], ((tpm_system->molecule+i)->atoms+j)->z);
			tpm_system->grid[idx][tpm_system->number_atom_in_grid[idx]]=(tpm_system->molecule+i)->atoms+j;
			tpm_system->number_atom_in_grid[idx]++;
		}
	}
/*	
	for(i=0; i<tpm_system->number_slice; i++)
	{
		for(j=0; j<tpm_system->number_slice; j++)
		{
			for(k=0; k<tpm_system->number_slice; k++)
			{
				printf("%d ",tpm_system->number_atom_in_grid[i*tpm_system->number_slice*tpm_system->number_slice+j*tpm_system->number_slice+k]);
			}
			printf("\n");
		}
		printf("\n");
	}
*/	
	return 1;
}

int	calculate_force(System *tpm_system)
{
	int i=0;
	int j=0;
	//Reset Accelerations
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			((tpm_system->molecule+i)->atoms+j)->ax = 0.0;
			((tpm_system->molecule+i)->atoms+j)->ay = 0.0;
			((tpm_system->molecule+i)->atoms+j)->az = 0.0;
		}
	}
	tpm_system->potential_energy = 0.0;
	CalculateBondForce(tpm_system);
	CalculateAngleForce(tpm_system);
	CalculatePairForce(tpm_system);
//	getchar();
	return 1;
}


int	integrate(System *tpm_system)
{
	int i=0;
	int j=0;
	double newx,newy,newz,dx,dy,dz;
	double dt=0.5;
	ATOM *atom;
//	double max = 0.0;
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			atom = (tpm_system->molecule+i)->atoms+j;
			newx = 2.0 *	atom->x - atom->oldx + atom->ax * dt * dt;
			newy = 2.0 *	atom->y - atom->oldy + atom->ay * dt * dt;
			newz = 2.0 *	atom->z - atom->oldz + atom->az * dt * dt;
			
			atom->oldx = atom->x;
			atom->oldy = atom->y;
			atom->oldz = atom->z;
			dx = newx - atom->oldx;
			dy = newy - atom->oldy;
			dz = newz - atom->oldz;
      get_minimun_image (&dx, &dy, &dz, tpm_system);

			atom->vx = dx/dt;
			atom->vy = dy/dt;
			atom->vz = dz/dt;
//			if(fabs(atom->ax) > 0.005) atom->ax=0.0;
//			if(fabs(atom->ay) > 0.005) atom->ay=0.0;
//			if(fabs(atom->az) > 0.005) atom->az=0.0;
										
			
			if(newx<0) newx+=tpm_system->dimension;
			if(newy<0) newy+=tpm_system->dimension;
			if(newz<0) newz+=tpm_system->dimension;
			if(newx>=tpm_system->dimension) newx-=tpm_system->dimension;
			if(newy>=tpm_system->dimension) newy-=tpm_system->dimension;
			if(newz>=tpm_system->dimension) newz-=tpm_system->dimension;
			atom->x = newx;
			atom->y = newy;
			atom->z = newz;
		}
	}
	return 1;
}

int	calculate_kinetic_energy(System *tpm_system)
{
	int i=0;
	int j=0;
	
	tpm_system->kinetic_energy = 0.0;
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			tpm_system->kinetic_energy += 	((tpm_system->molecule+i)->atoms+j)->vx*
																((tpm_system->molecule+i)->atoms+j)->vx*
																((tpm_system->molecule+i)->atoms+j)->mass*0.5;
			tpm_system->kinetic_energy += 	((tpm_system->molecule+i)->atoms+j)->vy*
																((tpm_system->molecule+i)->atoms+j)->vy*
																((tpm_system->molecule+i)->atoms+j)->mass*0.5;
			tpm_system->kinetic_energy += 	((tpm_system->molecule+i)->atoms+j)->vz*
																((tpm_system->molecule+i)->atoms+j)->vz*
																((tpm_system->molecule+i)->atoms+j)->mass*0.5;
		}
	}
	return 1;
}

int	update_system (System *tpm_system)
{
	integrate(tpm_system);
	create_grid_list(tpm_system);
	calculate_force(tpm_system);
//	RenderMolecules(tpm_system);
	release_grid_list(tpm_system);
	calculate_kinetic_energy(tpm_system);
	return 1;
}

void
integrate_relaxation (System *tpm_system)
{
	int i,j;
	double newx,newy,newz,dx,dy,dz;
	double dt=0.5;
	ATOM *atom;

	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			atom = (tpm_system->molecule+i)->atoms+j;
			
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
			
			get_minimun_image(&dx, &dy, &dz, tpm_system);

			atom->vx = dx/dt;
			atom->vy = dy/dt;
			atom->vz = dz/dt;
			
			if(newx<0) newx+=tpm_system->dimension;
			if(newy<0) newy+=tpm_system->dimension;
			if(newz<0) newz+=tpm_system->dimension;
			if(newx>=tpm_system->dimension) newx-=tpm_system->dimension;
			if(newy>=tpm_system->dimension) newy-=tpm_system->dimension;
			if(newz>=tpm_system->dimension) newz-=tpm_system->dimension;
			atom->x = newx;
			atom->y = newy;
			atom->z = newz;
		}
	}
}


int reset_velosity(System *tpm_system)
{
	int i,j;
	for(i=0; i<tpm_system->number_molecule; i++)
	{
		for(j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
			((tpm_system->molecule+i)->atoms+j)->vx=0.0;
			((tpm_system->molecule+i)->atoms+j)->vy=0.0;
			((tpm_system->molecule+i)->atoms+j)->vz=0.0;
			((tpm_system->molecule+i)->atoms+j)->ax=0.0;
			((tpm_system->molecule+i)->atoms+j)->ay=0.0;
			((tpm_system->molecule+i)->atoms+j)->az=0.0;
			((tpm_system->molecule+i)->atoms+j)->oldx=((tpm_system->molecule+i)->atoms+j)->x;
			((tpm_system->molecule+i)->atoms+j)->oldy=((tpm_system->molecule+i)->atoms+j)->y;
			((tpm_system->molecule+i)->atoms+j)->oldz=((tpm_system->molecule+i)->atoms+j)->z;
		}
	}
	return 1;
}



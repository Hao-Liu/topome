/*
 * /__ __X  _ \/  __\/  _ \/ \__/|/  __/ Copyright (C) 2010 - Liu Hao
 *   / \ | / \||  \/|| / \|| |\/|||  \*  
 *   | | | \_/||  __/| \_/|| |  |||  /_  force.c
 *   \_/ \____/\_/   \____/\_/  \|\____\
 * 
 * calculate forces and energies
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
#include "types.h"
void 
get_minimun_image (double *dx, double *dy, double *dz, System *tpm_system)
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


void
calculate_bond_force (System *tpm_system)
{
	int i, j;
	double dx, dy, dz, r, dr, force_coefficent, K, r0;
	for (i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule + i;
		for (j=0; j<mol->nBonds; j++)
		{
		  BOND *bond = (tpm_system->molecule+i)->bonds+j;
		  ATOM *atom1 = bond->atom1;
		  ATOM *atom2 = bond->atom2;
/*
// Used when manually input forcefield
			K  =  (tpm_system->bond_type+bond->type)->K;
			r0 =  (tpm_system->bond_type+bond->type)->r0;
*/

// Used When automatic generate forcefield 
			K  =  bond->k2;
			r0 =  bond->r0;

			dx =  atom1->x - atom2->x;
			dy =  atom1->y - atom2->y;
			dz =  atom1->z - atom2->z;
      
      get_minimun_image (&dx, &dy, &dz, tpm_system);

			r = sqrt(dx*dx+dy*dy+dz*dz);
			dr = r -r0;
			force_coefficent = 2.0 * K * dr / r *4.184e-4;
			
			tpm_system->potential_energy += K * dr *dr *4.184e-4;
			
			atom1->ax -= force_coefficent * dx / atom1->mass;
			atom1->ay -= force_coefficent * dy / atom1->mass;
			atom1->az -= force_coefficent * dz / atom1->mass;
			atom2->ax += force_coefficent * dx / atom2->mass;
			atom2->ay += force_coefficent * dy / atom2->mass;
			atom2->az += force_coefficent * dz / atom2->mass;
		}
	}
}

void
calculate_angle_force (System *tpm_system)
{
	int i,j;
	double K, theta0, dx1, dy1, dz1, rsq1, r1, dx2 ,dy2, dz2, rsq2, r2;
	double cos_theta, sin_theta, dtheta, tk, a, a11, a12, a22;
	double f1x, f1y, f1z, f3x, f3y, f3z;
	
	for (i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule+i;
		for (j=0; j<mol->nAngles; j++)
		{
		  ANGLE *angle = mol->angles + j;
/*
//Manual
			K =       (tpm_system->angle_type+angle->type)->K;
			// /180.0*3.14159265357989
			theta0 =  (tpm_system->angle_type+angle->type)->theta0 * 1.745329252e-2;
*/

//Automatic
			K =       angle->k2;
			// /180.0*3.14159265357989
			theta0 =  angle->theta0 * 1.745329252e-2;

			dx1 =  angle->atom1->x - angle->atom2->x;
			dy1 =  angle->atom1->y - angle->atom2->y;
			dz1 =  angle->atom1->z - angle->atom2->z;
      get_minimun_image (&dx1, &dy1, &dz1, tpm_system);
			rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
			r1 = sqrt (rsq1);

			dx2 =  angle->atom3->x - angle->atom2->x;
			dy2 =  angle->atom3->y - angle->atom2->y;
			dz2 =  angle->atom3->z - angle->atom2->z;
      get_minimun_image (&dx2, &dy2, &dz2, tpm_system);
			rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
			r2 = sqrt (rsq2);

			cos_theta = dx1*dx2 + dy1*dy2 + dz1*dz2;
			cos_theta /= r1*r2;
			if (cos_theta >  1.0) cos_theta =  1.0;
			if (cos_theta < -1.0) cos_theta = -1.0;
			
			sin_theta = sqrt (1.0 - cos_theta*cos_theta);
			if (sin_theta < 0.001) sin_theta = 0.001;
			dtheta = acos (cos_theta) - theta0;
			
			tk = K * dtheta;
			
			a   = -2.0 * tk / sin_theta * 4.184e-4;
			a11 = a*cos_theta / rsq1;
			a12 = -a / (r1*r2);
			a22 = a*cos_theta / rsq2;
			
			f1x = a11*dx1 + a12*dx2;
			f1y = a11*dy1 + a12*dy2;
			f1z = a11*dz1 + a12*dz2;
			f3x = a22*dx2 + a12*dx1;
			f3y = a22*dy2 + a12*dy1;
			f3z = a22*dz2 + a12*dz1;
			
			tpm_system->potential_energy += K * dtheta *dtheta * 4.184e-4;
			angle->atom1->ax += f1x / angle->atom1->mass;
			angle->atom1->ay += f1y / angle->atom1->mass;
			angle->atom1->az += f1z / angle->atom1->mass;
			angle->atom3->ax += f3x / angle->atom3->mass;
			angle->atom3->ay += f3y / angle->atom3->mass;
			angle->atom3->az += f3z / angle->atom3->mass;
			angle->atom2->ax -= (f1x+f3x) / angle->atom2->mass;
			angle->atom2->ay -= (f1y+f3y) / angle->atom2->mass;
			angle->atom2->az -= (f1z+f3z) / angle->atom2->mass;
		}
	}
}

//Pair force implemented grid scheme
void
calculate_pair_force_grid (System *tpm_system)
{
	int i,j,k,l,m,n,o,p,im,in,io;
	double dx, dy, dz, r2, invr2, invr, invr6;
	double lj_coefficient, coulomb_coefficient, force_coefficient;
	double coulomb_epsilon = 0.01;
	double lj_epsilon = 0.001;
	double soft = 0.001;

  double radius_cut_square = tpm_system->radius_cut*tpm_system->radius_cut;
  int number_slice_square = tpm_system->number_slice*tpm_system->number_slice;
  
	for (i=0; i<tpm_system->number_slice; i++)
	{
		for (j=0; j<tpm_system->number_slice; j++)
		{
			for (k=0; k<tpm_system->number_slice; k++)
			{
			  int grid_index1 = i * number_slice_square + 
			                    j * tpm_system->number_slice + k;

				for (l=0; l<tpm_system->number_atom_in_grid[grid_index1]; l++)
				{
					ATOM *atom1 = tpm_system->grid[grid_index1][l];
					for (im=i-1; im<i+2; im++)
					{
						for (in=j-1; in<j+2; in++)
						{
							for (io=k-1; io<k+2; io++)
							{
								m=im;
								n=in;
								o=io;
								if (m<0) m+=tpm_system->number_slice;
								if (n<0) n+=tpm_system->number_slice;
								if (o<0) o+=tpm_system->number_slice;
								if (m>=tpm_system->number_slice) m-=tpm_system->number_slice;
								if (n>=tpm_system->number_slice) n-=tpm_system->number_slice;
								if (o>=tpm_system->number_slice) o-=tpm_system->number_slice;
								
								int grid_index2 = m * number_slice_square +
								                  n * tpm_system->number_slice + o;
								                  
								for (p=0; p<tpm_system->number_atom_in_grid[grid_index2]; p++)
								{
									ATOM *atom2 = tpm_system->grid[grid_index2][p];
									if (atom1->mol == atom2->mol) continue;
									dx = atom2->x - atom1->x;
									dy = atom2->y - atom1->y;
									dz = atom2->z - atom1->z;
									
                  get_minimun_image (&dx, &dy, &dz, tpm_system);

									r2 = dx*dx + dy*dy + dz*dz + soft;
									if (r2 > radius_cut_square) continue; 
									invr2 = 1.0 / r2; 
									invr = sqrt (invr2);

									coulomb_coefficient = coulomb_epsilon * atom2->charge * 
									                      atom1->charge * invr * invr2;

									invr6 = invr2 * invr2 * invr2;
									
									lj_coefficient = 4.0 * lj_epsilon * (6.0 - 12.0 * invr6) *
									                 invr6 * invr2;
						
									force_coefficient = (lj_coefficient - coulomb_coefficient) / 
									                    atom1->mass;
									atom1->ax += force_coefficient * dx;
									atom1->ay += force_coefficient * dy;
									atom1->az += force_coefficient * dz;
									
									tpm_system->potential_energy -= 4.0 * lj_epsilon * 
									                                (1.0 - 1.0 * invr6) * invr6;
									tpm_system->potential_energy += coulomb_coefficient * 
									                                atom2->charge * 
									                                atom1->charge * invr;
								}
							}
						}
					}
				}
			}
		}
	}
}

//Global pair force 
void
calculate_pair_force_all (System *tpm_system)
{
	int i, j, k, l;
	double dx, dy, dz, r2, invr2, invr, invr6;
	double lj_coefficient, coulomb_coefficient, force_coefficient;
	double coulomb_epsilon = 0.01;
	double lj_epsilon = 0.001;
	double soft = 0.001;
  double radius_cut_square = tpm_system->radius_cut*tpm_system->radius_cut;

	for (i=0; i<tpm_system->number_molecule; i++)
	{
		for (j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
		  ATOM *atom1 = (tpm_system->molecule+i)->atoms + j;
			for (k=0; k<tpm_system->number_molecule; k++)
			{
				if (i == k) continue;
				
				for (l=0; l<(tpm_system->molecule+k)->nAtoms; l++)
				{
    		  ATOM *atom2 = (tpm_system->molecule+k)->atoms + j;
    		  
					dx = 	atom2->x - atom1->x;
					dy = 	atom2->y - atom1->y;
					dz = 	atom2->z - atom1->z;
          get_minimun_image (&dx, &dy, &dz, tpm_system);

					r2 = dx*dx + dy*dy + dz*dz + soft;
					if (r2 > radius_cut_square) continue; 
					invr2 = 1.0 / r2; 
					invr = sqrt (invr2);

					coulomb_coefficient = coulomb_epsilon * atom2->charge * 
					                      atom1->charge * invr * invr2;

					invr6 = invr2 * invr2 * invr2;
					
					lj_coefficient = 4.0 * lj_epsilon * (6.0 - 12.0 * invr6) *
					                 invr6 * invr2;
		
					force_coefficient = (lj_coefficient - coulomb_coefficient) / 
					                    atom1->mass;
					atom1->ax += force_coefficient * dx;
					atom1->ay += force_coefficient * dy;
					atom1->az += force_coefficient * dz;
					
					tpm_system->potential_energy -= 4.0 * lj_epsilon * 
					                                (1.0 - 1.0 * invr6) * invr6;
					tpm_system->potential_energy += coulomb_coefficient * 
					                                atom2->charge * 
					                                atom1->charge * invr;
				}
			}
		}
	}
}

void
release_grid_list (System *tpm_system)
{
	int i;
	for (i=0; i<tpm_system->number_grid; i++)
	{
		free (tpm_system->grid[i]);
		tpm_system->grid[i] = NULL;
	}
	free (tpm_system->grid);
  tpm_system->grid = NULL;
	free (tpm_system->number_atom_in_grid);
	tpm_system->number_atom_in_grid = NULL;
}

void
create_grid_list (System *tpm_system)
{
	int i,j,k,ix,iy,iz;
	
	tpm_system->number_slice = (int)floor (tpm_system->dimension / 
	                                       tpm_system->radius_cut);
	tpm_system->number_grid = tpm_system->number_slice *
	                          tpm_system->number_slice *
	                          tpm_system->number_slice;
	                          
	double inv_grid_dim = (double) (tpm_system->number_slice)
	                      /tpm_system->dimension;
	                      
	tpm_system->grid = malloc (tpm_system->number_grid * sizeof (ATOM**));
	tpm_system->number_atom_in_grid = malloc (tpm_system->number_grid *
	                                          sizeof (int));
	
	//get grid counts
	memset (tpm_system->number_atom_in_grid, 0, 
	        tpm_system->number_grid * sizeof (int));
	for (i=0; i<tpm_system->number_molecule; i++)
	{
		for (j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
		  ATOM *atom = (tpm_system->molecule+i)->atoms+j;
			ix = (int)floor (atom->x * inv_grid_dim);
			iy = (int)floor (atom->y * inv_grid_dim);
			iz = (int)floor (atom->z * inv_grid_dim);
			if (ix < 0) ix = 0;
			if (iy < 0) iy = 0;
			if (iz < 0) iz = 0;
			if (ix >= tpm_system->number_slice) ix = tpm_system->number_slice - 1;
			if (iy >= tpm_system->number_slice) iy = tpm_system->number_slice - 1;
			if (iz >= tpm_system->number_slice) iz = tpm_system->number_slice - 1;
			
			int idx = iz * tpm_system->number_slice * tpm_system->number_slice +
			            iy * tpm_system->number_slice + ix;
			tpm_system->number_atom_in_grid[idx]++;
		}
	}
	
	//allocate grid space
	for (i=0; i<tpm_system->number_slice; i++)
	{
		for (j=0; j<tpm_system->number_slice; j++)
		{
			for (k=0; k<tpm_system->number_slice; k++)
			{
				int idx = i * tpm_system->number_slice * tpm_system->number_slice +
				          j * tpm_system->number_slice + k;
				tpm_system->grid[idx] = malloc (tpm_system->number_atom_in_grid[idx] *
				                                sizeof (ATOM**));
			}
		}
	}

	//set grid
	memset (tpm_system->number_atom_in_grid, 0, 
	        tpm_system->number_grid * sizeof (int));
	for (i=0; i<tpm_system->number_molecule; i++)
	{
		for (j=0; j<(tpm_system->molecule+i)->nAtoms; j++)
		{
		  ATOM *atom = (tpm_system->molecule+i)->atoms+j;

			ix = (int)floor(atom->x * inv_grid_dim);
			iy = (int)floor(atom->y * inv_grid_dim);
			iz = (int)floor(atom->z * inv_grid_dim);
			if (ix < 0) ix = 0;
			if (iy < 0) iy = 0;
			if (iz < 0) iz = 0;
			if (ix >= tpm_system->number_slice) ix = tpm_system->number_slice - 1;
			if (iy >= tpm_system->number_slice) iy = tpm_system->number_slice - 1;
			if (iz >= tpm_system->number_slice) iz = tpm_system->number_slice - 1;
			int idx = iz * tpm_system->number_slice * tpm_system->number_slice +
			          iy * tpm_system->number_slice + ix;
			tpm_system->grid[idx][tpm_system->number_atom_in_grid[idx]] = 
			  (tpm_system->molecule+i)->atoms+j;
			tpm_system->number_atom_in_grid[idx]++;
		}
	}
}

void
calculate_aom_force (System *tpm_system)
{
  return;
}
void
calculate_force (System *tpm_system)
{
	int i, j;

	//reset accelerations
	for (i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule + i;
		for (j=0; j<mol->nAtoms; j++)
		{
			ATOM *atom = mol->atoms + j;
			atom->ax = 0.0;
			atom->ay = 0.0;
			atom->az = 0.0;
		}
	}
  
  //reset potential energy
	tpm_system->potential_energy = 0.0;

	calculate_bond_force (tpm_system);
	calculate_angle_force (tpm_system);
	calculate_pair_force_grid (tpm_system);
}


void
integrate (int relax, System *tpm_system)
{
	int i, j;
	double newx, newy, newz, dx, dy, dz;
	double dt     = tpm_system->time_interval;
	double dt2    = dt * dt;
	double inv_dt = 1.0 / dt;
  
	for (i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule+i;
		for (j=0; j<mol->nAtoms; j++)
		{
			ATOM *atom = mol->atoms+j;
			newx = 2.0*atom->x - atom->oldx + atom->ax*dt2;
			newy = 2.0*atom->y - atom->oldy + atom->ay*dt2;
			newz = 2.0*atom->z - atom->oldz + atom->az*dt2;
			
			atom->oldx = atom->x;
			atom->oldy = atom->y;
			atom->oldz = atom->z;
			
			if(relax)
			{
			  dx = newx - atom->oldx * tpm_system->relax_ratio;
			  dy = newy - atom->oldy * tpm_system->relax_ratio;
			  dz = newz - atom->oldz * tpm_system->relax_ratio;
			}
			else
			{
			  dx = newx - atom->oldx;
			  dy = newy - atom->oldy;
			  dz = newz - atom->oldz;
			}

      get_minimun_image (&dx, &dy, &dz, tpm_system);

			atom->vx = dx*inv_dt;
			atom->vy = dy*inv_dt;
			atom->vz = dz*inv_dt;
			
			if (newx < 0.0)                    newx += tpm_system->dimension;
			if (newy < 0.0)                    newy += tpm_system->dimension;
			if (newz < 0.0)                    newz += tpm_system->dimension;
			if (newx >= tpm_system->dimension) newx -= tpm_system->dimension;
			if (newy >= tpm_system->dimension) newy -= tpm_system->dimension;
			if (newz >= tpm_system->dimension) newz -= tpm_system->dimension;
			
			atom->x = newx;
			atom->y = newy;
			atom->z = newz;
		}
	}
}

void
calculate_kinetic_energy (System *tpm_system)
{
	int i, j;

	tpm_system->kinetic_energy = 0.0;

	for (i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule+i;
		for (j=0; j<mol->nAtoms; j++)
		{
		  ATOM *atom = mol->atoms+j;
			tpm_system->kinetic_energy += (atom->vx * atom->vx +
			                               atom->vy * atom->vy +
			                               atom->vz * atom->vz) 
			                                 * atom->mass * 0.5;
		}
	}
}

/**************Time Checker*************************
struct timeval tv, tickCount;

gettimeofday( &tv, NULL );
tickCount.tv_usec = tv.tv_usec;

gettimeofday( &tv, NULL );
tickCount.tv_usec = tv.tv_usec -tickCount.tv_usec;
printf("%ld ", tickCount.tv_usec);
tickCount.tv_usec = tv.tv_usec;
***************************************************/



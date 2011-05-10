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
get_minimum_image (double *dx, double *dy, double *dz, System *tpm_system)
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

			K  =  bond->k2;
			r0 =  bond->r0;

			dx =  atom1->x - atom2->x;
			dy =  atom1->y - atom2->y;
			dz =  atom1->z - atom2->z;
      
      get_minimum_image (&dx, &dy, &dz, tpm_system);

			r = sqrt(dx*dx+dy*dy+dz*dz);
			dr = r -r0;
			force_coefficent = 2.0 * K * dr / r * 4.184E-4;
			
			tpm_system->potential_energy += K * dr *dr;
			
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
			
			K =       angle->k2;
			// /180.0*3.14159265357989
			theta0 =  angle->theta0 * 1.745329252e-2;

			dx1 =  angle->atom1->x - angle->atom2->x;
			dy1 =  angle->atom1->y - angle->atom2->y;
			dz1 =  angle->atom1->z - angle->atom2->z;
      get_minimum_image (&dx1, &dy1, &dz1, tpm_system);
			rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
			r1 = sqrt (rsq1);

			dx2 =  angle->atom3->x - angle->atom2->x;
			dy2 =  angle->atom3->y - angle->atom2->y;
			dz2 =  angle->atom3->z - angle->atom2->z;
      get_minimum_image (&dx2, &dy2, &dz2, tpm_system);
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

void
calculate_dihedral_force (System *tpm_system)
{
	int i, j, k, n;
	double dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx2m,dy2m,dz2m;
	double f1x,f1y,f1z,f2x,f2y,f2z,f3x,f3y,f3z,f4x,f4y,f4z;
	double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
	double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
	double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
	double c,s,p,sx2,sy2,sz2;
	double K,phi0,cs;

	for (i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule+i;
		for (j=0; j<mol->nDihedrals; j++)
		{
		  DIHEDRAL *dihedral = mol->dihedrals + j;
			
			K =       dihedral->kphi;
			n =       dihedral->n;
			// /180.0*3.14159265357989
			phi0 =  dihedral->phi0 * 1.745329252e-2;

			dx1 =  dihedral->atom1->x - dihedral->atom2->x;
			dy1 =  dihedral->atom1->y - dihedral->atom2->y;
			dz1 =  dihedral->atom1->z - dihedral->atom2->z;
      get_minimum_image (&dx1, &dy1, &dz1, tpm_system);

			dx2 =  dihedral->atom3->x - dihedral->atom2->x;
			dy2 =  dihedral->atom3->y - dihedral->atom2->y;
			dz2 =  dihedral->atom3->z - dihedral->atom2->z;
      get_minimum_image (&dx2, &dy2, &dz2, tpm_system);

			dx2m = -dx2;
			dy2m = -dy2;
			dz2m = -dz2;
      get_minimum_image (&dx2m, &dy2m, &dz2m, tpm_system);
			
			dx3 =  dihedral->atom4->x - dihedral->atom3->x;
			dy3 =  dihedral->atom4->y - dihedral->atom3->y;
			dz3 =  dihedral->atom4->z - dihedral->atom3->z;
      get_minimum_image (&dx3, &dy3, &dz3, tpm_system);

			ax = dy1*dz2m - dz1*dy2m;
			ay = dz1*dx2m - dx1*dz2m;
			az = dx1*dy2m - dy1*dx2m;
			bx = dy3*dz2m - dz3*dy2m;
			by = dz3*dx2m - dx3*dz2m;
			bz = dx3*dy2m - dy3*dx2m;

			rasq = ax*ax + ay*ay + az*az;
			rbsq = bx*bx + by*by + bz*bz;
			rgsq = dx2m*dx2m + dy2m*dy2m + dz2m*dz2m;
			rg   = sqrt(rgsq);

			rginv = ra2inv = rb2inv = 0.0;
			if (rg > 0) rginv = 1.0/rg;
			if (rasq > 0) ra2inv = 1.0/rasq;
			if (rbsq > 0) rb2inv = 1.0/rbsq;
			rabinv = sqrt(ra2inv*rb2inv);

			c = (ax*bx + ay*by + az*bz) * rabinv;
			s = rg * rabinv * (ax*dx3 + ay*dy3 + az*dz3);

			//TODO: error check

			if (c > 1.0) c=1.0;
			if (c < -1.0) c=-1.0;

			p = 1.0;
			df1 = 0.0;

			for (k=0; k<n; k++)
			{
				ddf1 = p*c - df1*s;
				df1 = p*s + df1*c;
				p = ddf1;
			}

			if (phi0 - 0.0 < 0.01)
				cs = 1;
			else cs = -1;

			p = p*cs;
			df1 = df1*cs;
			df1 *= -n;
			p += 1.0;

			if (n == 0)
			{
				p = 1.0 + cs;
				df1 =0.0;
			}

			tpm_system->potential_energy += K * p * 4.184e-4;

			fg = dx1*dx2m + dy1*dy2m + dz1*dz2m;
			hg = dx3*dx2m + dy3*dy2m + dz3*dz2m;
			fga = fg*ra2inv*rginv;
			hgb = hg*rb2inv*rginv;
			gaa = -ra2inv*rg;
			gbb = rb2inv*rg;

			dtfx = gaa*ax;
			dtfy = gaa*ay;
			dtfz = gaa*az;
			dtgx = fga*ax -hgb*bx;
			dtgy = fga*ay -hgb*by;
			dtgz = fga*az -hgb*bz;
			dthx = gbb*bx;
			dthy = gbb*by;
			dthz = gbb*bz;

			df = -K * df1 * 4.184e-4;

			sx2 = df*dtgx;
			sy2 = df*dtgy;
			sz2 = df*dtgz;
			
			f1x = df*dtfx;
			f1y = df*dtfy;
			f1z = df*dtfz;

			f2x = sx2 -f1x;
			f2y = sy2 -f1y;
			f2z = sz2 -f1z;

			f4x = df*dthx;
			f4y = df*dthy;
			f4z = df*dthz;

			f3x = -sx2 -f4x;
			f3y = -sy2 -f4y;
			f3z = -sz2 -f4z;

			dihedral->atom1->ax += f1x / dihedral->atom1->mass;
			dihedral->atom1->ay += f1y / dihedral->atom1->mass;
			dihedral->atom1->az += f1z / dihedral->atom1->mass;
			dihedral->atom3->ax += f3x / dihedral->atom3->mass;
			dihedral->atom3->ay += f3y / dihedral->atom3->mass;
			dihedral->atom3->az += f3z / dihedral->atom3->mass;
			dihedral->atom2->ax += f2x / dihedral->atom2->mass;
			dihedral->atom2->ay += f2y / dihedral->atom2->mass;
			dihedral->atom2->az += f2z / dihedral->atom2->mass;
			dihedral->atom4->ax += f4x / dihedral->atom4->mass;
			dihedral->atom4->ay += f4y / dihedral->atom4->mass;
			dihedral->atom4->az += f4z / dihedral->atom4->mass;
		}
	}
}

void
calculate_improper_force (System *tpm_system)
{
  int i,j,n;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double f1[3],f2[3],f3[3],f4[3];
  double sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
  double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
  double c2mag,sc1,sc2,s1,s2,s12,c,p,pd,rc2,a,a11,a22;
  double a33,a12,a13,a23,sx2,sy2,sz2;
	double K, chi0;


	for (i=0; i<tpm_system->number_molecule; i++)
	{
	  MOLECULE *mol = tpm_system->molecule+i;
		for (j=0; j<mol->nImpropers; j++)
		{
		  IMPROPER *improper = mol->impropers + j;
			
			K =       improper->kchi;
			n =       improper->n;
			// /180.0*3.14159265357989
			chi0 =  improper->chi0 * 1.745329252e-2;

			// 1st bond

			vb1x =  improper->atom1->x - improper->atom2->x;
			vb1y =  improper->atom1->y - improper->atom2->y;
			vb1z =  improper->atom1->z - improper->atom2->z;
      get_minimum_image (&vb1x, &vb1y, &vb1z, tpm_system);

			// 2nd bond

			vb2x =  improper->atom3->x - improper->atom2->x;
			vb2y =  improper->atom3->y - improper->atom2->y;
			vb2z =  improper->atom3->z - improper->atom2->z;
      get_minimum_image (&vb1x, &vb1y, &vb1z, tpm_system);

			vb2xm = -vb2x;
			vb2ym = -vb2y;
			vb2zm = -vb2z;
      get_minimum_image (&vb1x, &vb1y, &vb1z, tpm_system);

			// 3rd bond

			vb3x =  improper->atom4->x - improper->atom3->x;
			vb3y =  improper->atom4->y - improper->atom3->y;
			vb3z =  improper->atom4->z - improper->atom3->z;
      get_minimum_image (&vb3x, &vb3y, &vb3z, tpm_system);

			// c0 calculation

			sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
			sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
			sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

			rb1 = sqrt(sb1);
			rb3 = sqrt(sb3);

			c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

			// 1st and 2nd angle

			b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
			b1mag = sqrt(b1mag2);
			b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
			b2mag = sqrt(b2mag2);
			b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
			b3mag = sqrt(b3mag2);

			ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
			r12c1 = 1.0 / (b1mag*b2mag);
			c1mag = ctmp * r12c1;

			ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
			r12c2 = 1.0 / (b2mag*b3mag);
			c2mag = ctmp * r12c2;

			// cos and sin of 2 angles and final c

			sc1 = sqrt(1.0 - c1mag*c1mag);
			if (sc1 < 0.001) sc1 = 0.001;
			sc1 = 1.0/sc1;

			sc2 = sqrt(1.0 - c2mag*c2mag);
			if (sc2 < 0.001) sc2 = 0.001;
			sc2 = 1.0/sc2;

			s1 = sc1 * sc1;
			s2 = sc2 * sc2;
			s12 = sc1 * sc2;
			c = (c0 + c1mag*c2mag) * s12;
			
			//TODO: error check
			if (c > 1.0) c = 1.0;
			if (c < -1.0) c = -1.0;

			// force & energy
			// p = 1 + cos(n*phi) for d = 1
			// p = 1 - cos(n*phi) for d = -1
			// pd = dp/dc / 2

			if (n == 2) {
				p = 2.0*c*c;
				pd = 2.0*c;
			} else if (n == 3) {
				rc2 = c*c;
				p = (4.0*rc2-3.0)*c + 1.0;
				pd = 6.0*rc2 - 1.5;
			} else if (n == 4) {
				rc2 = c*c;
				p = 8.0*(rc2-1)*rc2 + 2.0;
				pd = (16.0*rc2-8.0)*c;
			} else if (n == 6) {
				rc2 = c*c;
				p = ((32.0*rc2-48.0)*rc2 + 18.0)*rc2;
				pd = (96.0*(rc2-1.0)*rc2 + 18.0)*c;
			} else if (n == 1) {
				p = c + 1.0;
				pd = 0.5;
			} else if (n == 5) {
				rc2 = c*c;
				p = ((16.0*rc2-20.0)*rc2 + 5.0)*c + 1.0;
				pd = (40.0*rc2-30.0)*rc2 + 2.5;
			} else if (n == 0) {
				p = 2.0;
				pd = 0.0;
			}

			if (chi0 - 0.0 > 0.01) {
				p = 2.0 - p;
				pd = -pd;
			}

			tpm_system->potential_energy += K * p * 4.184e-4;

			a = 2.0 * K * pd  * 4.184E-4;
			c = c * a;
			s12 = s12 * a;
			a11 = c*sb1*s1;
			a22 = -sb2*(2.0*c0*s12 - c*(s1+s2));
			a33 = c*sb3*s2;
			a12 = -r12c1*(c1mag*c*s1 + c2mag*s12);
			a13 = -rb1*rb3*s12;
			a23 = r12c2*(c2mag*c*s2 + c1mag*s12);

			sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
			sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
			sz2  = a12*vb1z + a22*vb2z + a23*vb3z;

			f1[0] = a11*vb1x + a12*vb2x + a13*vb3x;
			f1[1] = a11*vb1y + a12*vb2y + a13*vb3y;
			f1[2] = a11*vb1z + a12*vb2z + a13*vb3z;

			f2[0] = -sx2 - f1[0];
			f2[1] = -sy2 - f1[1];
			f2[2] = -sz2 - f1[2];

			f4[0] = a13*vb1x + a23*vb2x + a33*vb3x;
			f4[1] = a13*vb1y + a23*vb2y + a33*vb3y;
			f4[2] = a13*vb1z + a23*vb2z + a33*vb3z;

			f3[0] = sx2 - f4[0];
			f3[1] = sy2 - f4[1];
			f3[2] = sz2 - f4[2];

			// apply force to each of 4 atoms

			improper->atom1->ax += f1[0] / improper->atom1->mass;
			improper->atom1->ay += f1[1] / improper->atom1->mass;
			improper->atom1->az += f1[2] / improper->atom1->mass;
			improper->atom3->ax += f3[0] / improper->atom3->mass;
			improper->atom3->ay += f3[1] / improper->atom3->mass;
			improper->atom3->az += f3[2] / improper->atom3->mass;
			improper->atom2->ax += f2[0] / improper->atom2->mass;
			improper->atom2->ay += f2[1] / improper->atom2->mass;
			improper->atom2->az += f2[2] / improper->atom2->mass;
			improper->atom4->ax += f4[0] / improper->atom4->mass;
			improper->atom4->ay += f4[1] / improper->atom4->mass;
			improper->atom4->az += f4[2] / improper->atom4->mass;
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
									
                  get_minimum_image (&dx, &dy, &dz, tpm_system);

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
          get_minimum_image (&dx, &dy, &dz, tpm_system);

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
	calculate_dihedral_force (tpm_system);
	calculate_improper_force (tpm_system);
	calculate_pair_force_grid (tpm_system);
}

void
relax_integrate (System *tpm_system)
{
	int i, j;
	for (i=0; i<tpm_system->number_molecule; i++)
	{
		MOLECULE *mol = tpm_system->molecule+i;
		for (j=0; j<mol->nAtoms; j++)
		{
			ATOM *atom = mol->atoms+j;
			atom->x = atom->x * tpm_system->relax_ratio;
			atom->y = atom->y * tpm_system->relax_ratio;
			atom->z = atom->z * tpm_system->relax_ratio;
			
			atom->oldx = atom->oldx * tpm_system->relax_ratio;
			atom->oldy = atom->oldy * tpm_system->relax_ratio;
			atom->oldz = atom->oldz * tpm_system->relax_ratio;
		}
	}
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
			
			if (relax)
			{
				atom->x = atom->x * tpm_system->relax_ratio;
				atom->y = atom->y * tpm_system->relax_ratio;
				atom->z = atom->z * tpm_system->relax_ratio;
				
				atom->oldx = atom->oldx * tpm_system->relax_ratio;
				atom->oldy = atom->oldy * tpm_system->relax_ratio;
				atom->oldz = atom->oldz * tpm_system->relax_ratio;
			}

			newx = 2.0*atom->x - atom->oldx + atom->ax*dt2;
			newy = 2.0*atom->y - atom->oldy + atom->ay*dt2;
			newz = 2.0*atom->z - atom->oldz + atom->az*dt2;
			
			atom->oldx = atom->x;
			atom->oldy = atom->y;
			atom->oldz = atom->z;
			  
			dx = newx - atom->oldx;
			dy = newy - atom->oldy;
			dz = newz - atom->oldz;

      get_minimum_image (&dx, &dy, &dz, tpm_system);

			atom->vx = dx*inv_dt;
			atom->vy = dy*inv_dt;
			atom->vz = dz*inv_dt;
			
			if(relax)
			{
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
			else 
			{
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



#include "forcefield.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void
get_define (FILE *fp, long pos, ForceField* ff)
{
  ff->num_force_type++;
  ff->force_type = realloc (ff->force_type, 
                            sizeof (ForceType)*ff->num_force_type);
  ForceType *ft = ff->force_type + ff->num_force_type - 1;
  memset (ft, 0, sizeof (ForceType));

  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);
  sscanf (buffer+1, "%s", ft->name);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ft->num_function++;
        ft->function = realloc (ft->function, 
                                sizeof (ForceFunction) * ft->num_function);
        ForceFunction *func = ft->function + ft->num_function - 1;
        memset (func, 0, sizeof (ForceFunction) );

        sscanf (buffer, "%s %d %s %s %s", func->version, 
                                          &func->reference, 
                                          func->function, 
                                          func->label1, 
                                          func->label2);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }

  fseek (fp, next_section, SEEK_SET);
}

void
get_atom_type (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_atom_type++;
        ff->atom_type = realloc (ff->atom_type, 
                                 sizeof (AtomType) * ff->num_atom_type);
        AtomType *at = ff->atom_type + ff->num_atom_type - 1;
        memset (at, 0, sizeof (AtomType) );

        sscanf (buffer, "%s %d %s %f %s %d", 
          at->version, 
          &at->reference, 
          at->type, 
          &at->mass, 
          at->element, 
          &at->connection);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_equivalence (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_equivalence++;
        ff->equivalence = realloc (ff->equivalence, 
                                   sizeof (Equivalence) * ff->num_equivalence);
        Equivalence *eq = ff->equivalence + ff->num_equivalence - 1;
        memset (eq, 0, sizeof (Equivalence) );

        sscanf (buffer, "%s %d %s %s %s %s %s %s", 
          eq->version, 
          &eq->reference, 
          eq->type, 
          eq->nonbond, 
          eq->bond, 
          eq->angle, 
          eq->torsion, 
          eq->out_of_plane);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_auto_equivalence (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_auto_equivalence++;
        ff->auto_equivalence = realloc (ff->auto_equivalence, 
                                   sizeof (AutoEquivalence) * ff->num_auto_equivalence);
        AutoEquivalence *aeq = ff->auto_equivalence + ff->num_auto_equivalence - 1;
        memset (aeq, 0, sizeof (AutoEquivalence) );

        sscanf (buffer, "%s %d %s %s %s %s %s %s %s %s %s %s", 
          aeq->version, 
          &aeq->reference, 
          aeq->type, 
          aeq->nonbond, 
          aeq->bond_increment, 
          aeq->bond, 
          aeq->angle_end_atom, 
          aeq->angle_apex_atom, 
          aeq->torsion_end_atom, 
          aeq->torsion_center_atom, 
          aeq->out_of_plane_end_atom,
          aeq->out_of_plane_center_atom);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_hbond_definition (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  char entry[20] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  ff->hbond_definition = malloc (sizeof (HBondDefinition));
  HBondDefinition *hb = ff->hbond_definition;
  memset (hb, 0, sizeof (HBondDefinition) );
  
  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        sscanf (buffer, "%s %d %s", 
          hb->version, 
          &hb->reference, 
          entry);
        if (!strcmp (entry, "distance"))
        {
          sscanf (buffer, "%s %d %s %f", 
              hb->version,
              &hb->reference,
              entry,
              &hb->distance);
          next_section = ftell (fp);
          continue;
        }
        if (!strcmp (entry, "angle"))
        {
          sscanf (buffer, "%s %d %s %f", 
              hb->version,
              &hb->reference,
              entry,
              &hb->angle);
          next_section = ftell (fp);
          continue;
        }
        if (!strcmp (entry, "donors"))
        {
          sscanf (buffer, "%s %d %s %s %s %s %s", //FIXME not generic for all frc 
              hb->version,
              &hb->reference,
              entry,
              hb->donor[0],
              hb->donor[1],
              hb->donor[2],
              hb->donor[3]);
          next_section = ftell (fp);
          continue;
        }
        if (!strcmp (entry, "acceptors"))
        {
          sscanf (buffer, "%s %d %s %s %s %s %s %s", //FIXME not generic for all frc 
              hb->version,
              &hb->reference,
              entry,
              hb->acceptor[0],
              hb->acceptor[1],
              hb->acceptor[2],
              hb->acceptor[3],
              hb->acceptor[4]);
          next_section = ftell (fp);
          continue;
        }
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_morse_bond (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_morse_bond++;
        ff->morse_bond = realloc (ff->morse_bond, 
                                 sizeof (MorseBond) * ff->num_morse_bond);
        MorseBond *mb = ff->morse_bond + ff->num_morse_bond - 1;
        memset (mb, 0, sizeof (MorseBond) );

        sscanf (buffer, "%s %d %s %s %f %f %f", 
          mb->version, 
          &mb->reference, 
          mb->i, 
          mb->j, 
          &mb->r0, 
          &mb->d,
          &mb->alpha);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_quadratic_bond (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_quadratic_bond++;
        ff->quadratic_bond = realloc (ff->quadratic_bond, 
                                 sizeof (QuadraticBond) * ff->num_quadratic_bond);
        QuadraticBond *qb = ff->quadratic_bond + ff->num_quadratic_bond - 1;
        memset (qb, 0, sizeof (QuadraticBond) );

        sscanf (buffer, "%s %d %s %s %f %f", 
          qb->version, 
          &qb->reference, 
          qb->i, 
          qb->j, 
          &qb->r0, 
          &qb->k2);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_quadratic_angle (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_quadratic_angle++;
        ff->quadratic_angle = realloc (ff->quadratic_angle, 
                                 sizeof (QuadraticAngle) * ff->num_quadratic_angle);
        QuadraticAngle *qa = ff->quadratic_angle + ff->num_quadratic_angle - 1;
        memset (qa, 0, sizeof (QuadraticAngle) );

        sscanf (buffer, "%s %d %s %s %s %f %f", 
          qa->version, 
          &qa->reference, 
          qa->i, 
          qa->j, 
          qa->k, 
          &qa->theta0, 
          &qa->k2);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_bond_bond (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_bond_bond++;
        ff->bond_bond = realloc (ff->bond_bond, 
                                 sizeof (BondBond) * ff->num_bond_bond);
        BondBond *bb = ff->bond_bond + ff->num_bond_bond - 1;
        memset (bb, 0, sizeof (BondBond) );

        sscanf (buffer, "%s %d %s %s %s %f", 
          bb->version, 
          &bb->reference, 
          bb->i, 
          bb->j, 
          bb->k, 
          &bb->kbb);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_bond_angle (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_bond_angle++;
        ff->bond_angle = realloc (ff->bond_angle, 
                                 sizeof (BondAngle) * ff->num_bond_angle);
        BondAngle *ba = ff->bond_angle + ff->num_bond_angle - 1;
        memset (ba, 0, sizeof (BondAngle) );

        sscanf (buffer, "%s %d %s %s %s %f %f", 
          ba->version, 
          &ba->reference, 
          ba->i, 
          ba->j, 
          ba->k, 
          &ba->kbt1,
          &ba->kbt2);
        if (!strcmp (ba->i, ba->k))
        {
          ba->kbt2 = ba->kbt1;
        }
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_torsion_1 (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_torsion_1++;
        ff->torsion_1 = realloc (ff->torsion_1, 
                                 sizeof (Torsion1) * ff->num_torsion_1);
        Torsion1 *t1 = ff->torsion_1 + ff->num_torsion_1 - 1;
        memset (t1, 0, sizeof (Torsion1) );

        sscanf (buffer, "%s %d %s %s %s %s %f %f %f", 
          t1->version, 
          &t1->reference, 
          t1->i, 
          t1->j, 
          t1->k, 
          t1->l, 
          &t1->kphi,
          &t1->n,
          &t1->phi0);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_angle_angle_torsion_1 (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_angle_angle_torsion_1++;
        ff->angle_angle_torsion_1 = realloc (ff->angle_angle_torsion_1, 
                                 sizeof (AngleAngleTorsion1) * ff->num_angle_angle_torsion_1);
        AngleAngleTorsion1 *aat1 = ff->angle_angle_torsion_1 + ff->num_angle_angle_torsion_1 - 1;
        memset (aat1, 0, sizeof (AngleAngleTorsion1) );

        sscanf (buffer, "%s %d %s %s %s %s %f", 
          aat1->version, 
          &aat1->reference, 
          aat1->i, 
          aat1->j, 
          aat1->k, 
          aat1->l, 
          &aat1->kaat);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_out_of_plane (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_out_of_plane++;
        ff->out_of_plane = realloc (ff->out_of_plane, 
                                 sizeof (OutOfPlane) * ff->num_out_of_plane);
        OutOfPlane *oop = ff->out_of_plane + ff->num_out_of_plane - 1;
        memset (oop, 0, sizeof (OutOfPlane) );

        sscanf (buffer, "%s %d %s %s %s %s %f %f %f", 
          oop->version, 
          &oop->reference, 
          oop->i, 
          oop->j, 
          oop->k, 
          oop->l, 
          &oop->kchi,
          &oop->n,
          &oop->chi0);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_out_of_plane_out_of_plane (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_out_of_plane_out_of_plane++;
        ff->out_of_plane_out_of_plane = realloc (ff->out_of_plane_out_of_plane, 
                                 sizeof (OutOfPlaneOutOfPlane) * ff->num_out_of_plane_out_of_plane);
        OutOfPlaneOutOfPlane *oopoop = ff->out_of_plane_out_of_plane + ff->num_out_of_plane_out_of_plane - 1;
        memset (oopoop, 0, sizeof (OutOfPlaneOutOfPlane) );

        sscanf (buffer, "%s %d %s %s %s %s %f", 
          oopoop->version, 
          &oopoop->reference, 
          oopoop->i, 
          oopoop->j, 
          oopoop->k, 
          oopoop->l, 
          &oopoop->koo);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_angle_angle (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_angle_angle++;
        ff->angle_angle = realloc (ff->angle_angle, 
                                 sizeof (AngleAngle) * ff->num_angle_angle);
        AngleAngle *aa = ff->angle_angle + ff->num_angle_angle - 1;
        memset (aa, 0, sizeof (AngleAngle) );

        sscanf (buffer, "%s %d %s %s %s %s %f", 
          aa->version, 
          &aa->reference, 
          aa->i, 
          aa->j, 
          aa->k, 
          aa->l, 
          &aa->kaa);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_nonbond_lj (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_nonbond_lj++;
        ff->nonbond_lj = realloc (ff->nonbond_lj, 
                                 sizeof (NonbondLJ) * ff->num_nonbond_lj);
        NonbondLJ *nblj = ff->nonbond_lj + ff->num_nonbond_lj - 1;
        memset (nblj, 0, sizeof (NonbondLJ) );

        sscanf (buffer, "%s %d %s %f %f", 
          nblj->version, 
          &nblj->reference, 
          nblj->i, 
          &nblj->a, 
          &nblj->b);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_bond_increments (FILE *fp, long pos, ForceField* ff)
{
  char buffer[2000] = {0};
  fseek (fp, pos, SEEK_SET);

  fgets (buffer, sizeof (buffer), fp);

  long next_section = ftell (fp);
  while (fgets (buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case ' ':
        ff->num_bond_increments++;
        ff->bond_increments = realloc (ff->bond_increments, 
                                 sizeof (BondIncrements) * ff->num_bond_increments);
        BondIncrements *bi = ff->bond_increments + ff->num_bond_increments - 1;
        memset (bi, 0, sizeof (BondIncrements) );

        sscanf (buffer, "%s %d %s %s %f %f", 
          bi->version, 
          &bi->reference, 
          bi->i, 
          bi->j, 
          &bi->dij, 
          &bi->dji);
        next_section = ftell (fp);
        continue;
      case '#':
        break;
      default:
        next_section = ftell (fp);
        continue;
    }
    break;
  }
  fseek (fp, next_section, SEEK_SET);
}

void
get_section (FILE *fp, long pos, ForceField* ff)
{
  fseek (fp, pos, SEEK_SET);
  
  char buffer[2000] = {0};
  char section_name[200] = {0};

  fgets (buffer, sizeof (buffer), fp);
  if (buffer[0] == '#')
  {
    sscanf (buffer+1, "%s", section_name);
    if ( !strcmp (section_name, "version"))
    {
      return;
    }
    if ( !strcmp (section_name, "define"))
    {
      get_define (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "atom_types"))
    {
      get_atom_type (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "equivalence"))
    {
      get_equivalence (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "auto_equivalence"))
    {
      get_auto_equivalence (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "hbond_definition"))
    {
      get_hbond_definition (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "morse_bond"))
    {
      get_morse_bond (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "quadratic_bond"))
    {
      get_quadratic_bond (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "quadratic_angle"))
    {
      get_quadratic_angle (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "bond-bond"))
    {
      get_bond_bond (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "bond-angle"))
    {
      get_bond_angle (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "torsion_1"))
    {
      get_torsion_1 (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "angle-angle-torsion_1"))
    {
      get_angle_angle_torsion_1 (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "out_of_plane"))
    {
      get_out_of_plane (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "out_of_plane-out_of_plane"))
    {
      get_out_of_plane_out_of_plane (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "angle-angle"))
    {
      get_angle_angle (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "nonbond(12-6)"))
    {
      get_nonbond_lj (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "bond_increments"))
    {
      get_bond_increments (fp, pos, ff);
      return;
    }
    if ( !strcmp (section_name, "reference"))
    {
      return;
    }
    if ( !strcmp (section_name, "end"))
    {
      return;
    }
    else
    {
      printf("No %s Defined!\n", section_name);
      return;
    }
  }
  else
  {
    return; //FIXME error handling;
  }
}

void
init_forcefield(char *filename, ForceField* forcefield)
{
  FILE *fp = fopen (filename, "r");
  if (!fp)
    return; //FIXME: error handling;
  
  char buffer[2000];
  
  long section = ftell (fp);
  while (fgets(buffer, sizeof (buffer), fp))
  {
    switch (buffer[0])
    {
      case '#':
        get_section (fp, section, forcefield);
      default:
        section = ftell (fp);
        continue;
    }
  }
  
  fclose (fp);
}

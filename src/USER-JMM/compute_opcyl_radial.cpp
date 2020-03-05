/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"
#include "compute_opcyl_radial.h"

using namespace LAMMPS_NS;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Compute an order parameter which is the standard deviation of 
// the radial positions of each atom. It is assumed that the central
// axis is defined by z = 0.
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ComputeOPCylRadial::ComputeOPCylRadial(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  //std::cout << "ComputeOPCylRadial::ComputeOPCylRadial() has been called." << std::endl;
  // Arguments [3]:
  // 0. ID
  // 1. group-ID
  // 2. opcyl/radial
  //

  if (narg != 3) error->all(FLERR,"Illegal compute opcyl/radial command");

  scalar_flag = 1;
  vector_flag=0;
  array_flag=0;
  extscalar = 0;
  tempflag=0;
  extarray=0;
  
}


/* ---------------------------------------------------------------------- */

ComputeOPCylRadial::~ComputeOPCylRadial()
{
}


/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

double ComputeOPCylRadial::compute_scalar()
{
  
  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double r;
  double sum[2] = {0.0,0.0}, sumall[2];
  int ntotal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      r      = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
      sum[0] += r;
      sum[1] += r*r;
    }
    else {
    }
  }

  MPI_Allreduce(&sum,&sumall,2,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&nlocal,&ntotal,1,MPI_INT,MPI_SUM,world);
  scalar = ntotal*sumall[1] - sumall[0]*sumall[0];
  if (scalar > 0.0) {
    scalar = sqrt( scalar ) / ntotal;
  }
  else {
    scalar = 0.0;
  }

  return scalar;
}


/* ---------------------------------------------------------------------- */


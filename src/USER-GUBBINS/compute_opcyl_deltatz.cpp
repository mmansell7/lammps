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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>
#include "compute_opcyl_deltatz.h"
#include "math_const.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
//#include "force.h"
//#include "pair.h"
#include "comm.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Compute an order parameter which is the standard deviation of 
// the difference in azimuthal angle divided by the z-dimension separations 
// of each adjacent atom pair. It is assumed that the central axis is 
// defined by z = 0.
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* ---------------------------------------------------------------------- */

ComputeOPCylDeltaTZ::ComputeOPCylDeltaTZ(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // Arguments [3]:
  // 0. ID
  // 1. group-ID
  // 2. opcyl/deltatz
  //

  if (narg != 3) error->all(FLERR,"Illegal compute opcyl/deltatz command");

  scalar_flag = 1;
  vector_flag=0;
  array_flag=0;
  extscalar = 0;
  tempflag=0;
  extarray=0;
  
}

/* ---------------------------------------------------------------------- */

ComputeOPCylDeltaTZ::~ComputeOPCylDeltaTZ()
{
}

/* ---------------------------------------------------------------------- */

void ComputeOPCylDeltaTZ::init()
{
  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"opcyl/deltatz") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute opcyl/deltatz");
  
}

/* ---------------------------------------------------------------------- */

void ComputeOPCylDeltaTZ::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeOPCylDeltaTZ::compute_scalar()
{
  int i;
  double zl,dz,dphi,s;
  double sum = 0, sumsq = 0;
  
  invoked_scalar = update->ntimestep;
  
  double **x     = atom->x;
  int *mask      = atom->mask;
  int nlocal     = atom->nlocal;
  bigint natoms  = atom->natoms;
 
  int *nlocals   = new int[comm->nprocs];
  int *displs    = new int[comm->nprocs];
  tzlocal        = new double[2*nlocal];
  tzglobal       = new double[2*natoms];
  
  // Give all processes the number of atoms local to all other processes
  nlocal = 2*nlocal;
  MPI_Allgather(&nlocal,1,MPI_INT,nlocals,1,MPI_INT,world);
  
  // Set up vector of displacements so that when z coordinates are gathered,
  //   each processor knows where in the global z vector the z coordinates
  //   from each processor will begin.
  displs[0] = 0;
  for ( i = 1; i < comm->nprocs; i++ ) {
    displs[i] = displs[i-1] + nlocals[i-1];
  }
  
  // Assign z and phi coordinates for local atoms to tzlocal, the gather all the tzlocals
  //   into a "shared" vector of all atoms' z and phi coordinates, called tzglobal.
  for ( i = 0; i < nlocal/2; i++ ) {
    tzlocal[2*i] = x[i][2];
    tzlocal[2*i+1] = atan2(x[i][1],x[i][0]);
  }
  MPI_Allgatherv(tzlocal,nlocal,MPI_DOUBLE,tzglobal,nlocals,displs,MPI_DOUBLE,world);
  
  // The root process sorts the z-positions in order, then uses the differences
  //   between adjacent indices in the vector to obtain a vector of delta-z
  //   values.
  if ( comm->me == 0 ) {
    sortdoubles(tzglobal,natoms);
    zl = domain->boxhi[2] - domain->boxlo[2];   // z-length of simulation cell
    dz = tzglobal[0] - tzglobal[2*natoms-2];            // Wrap the first
    dz = dz - zl * floor( dz / zl); //  delta-z value.
    dphi = tzglobal[1] - tzglobal[2*natoms-1];
    dphi = remainder(dphi,MY_2PI);
    s = fabs(dphi/dz);
    sum = s;
    sumsq = s*s;
    for ( i = 1; i < natoms; i++ ) {
      dz = tzglobal[2*i] - tzglobal[2*i-2];
      dphi = tzglobal[2*i+1] - tzglobal[2*i-1];
      dphi = remainder(dphi,MY_2PI);
      s = fabs(dphi/dz);
      sum += s;
      sumsq += s*s;
    }
    scalar = natoms*sumsq - sum*sum;            // A preliminary value
    // Ensure small negative values that are essentially equal to 0 do not
    //   lead to nan result from the sqrt.
    if ( scalar < 0.0 && scalar > -1.0E-8 ) scalar = 0.0;
    scalar = sqrt( scalar ) / natoms; // Get standard dev from the preliminary value
  }
  
  // Broadcast the value of scalar calculated by the root process to the other processes
  MPI_Bcast(&scalar,1,MPI_DOUBLE,0,world);
  
  // Be sure to clean up memory to avoid memory leaks
  delete[] tzlocal;
  delete[] tzglobal;
  delete[] nlocals;
  delete[] displs;
  
  return scalar;
}


/* ---------------------------------------------------------------------- */
// A function that sorts a vector (tz) containing n doubles
void ComputeOPCylDeltaTZ::sortdoubles(double *tz,bigint n) {
  int i,j;
  double a;
  
  for (i = 0; i < 2*n-2; i += 2) {
    for (j = i + 2; j < 2*n; j += 2) {
      if (tz[i] > tz[j]) {
        a = tz[i];
        tz[i] = tz[j];
        tz[j] = a;
        a = tz[i+1];
        tz[i+1] = tz[j+1];
        tz[j+1] = a;
      }
    }
  }
}



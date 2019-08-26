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
#include "compute_opcyl_deltaz.h"
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

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Compute an order parameter which is the standard deviation of 
// the z-dimension separations of each adjacent atom pair. It is
// assumed that the central axis is defined by z = 0.
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* ---------------------------------------------------------------------- */

ComputeOPCylDeltaZ::ComputeOPCylDeltaZ(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // Arguments [3]:
  // 0. ID
  // 1. group-ID
  // 2. opcyl/deltaz
  //

  if (narg != 3) error->all(FLERR,"Illegal compute opcyl/deltaz command");

  scalar_flag = 1;
  vector_flag=0;
  array_flag=0;
  extscalar = 0;
  tempflag=0;
  extarray=0;
  
}

/* ---------------------------------------------------------------------- */

ComputeOPCylDeltaZ::~ComputeOPCylDeltaZ()
{
}

/* ---------------------------------------------------------------------- */

void ComputeOPCylDeltaZ::init()
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
    if (strcmp(modify->compute[i]->style,"opcyl/deltaz") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute opcyl/deltaz");
  
}

/* ---------------------------------------------------------------------- */

void ComputeOPCylDeltaZ::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeOPCylDeltaZ::compute_scalar()
{
  int i;
  double zl;
  double sum = 0, sumsq = 0;
  
  invoked_scalar = update->ntimestep;
  
  double **x     = atom->x;
  int *mask      = atom->mask;
  int nlocal     = atom->nlocal;
  bigint natoms  = atom->natoms;
 
  int *nlocals   = new int[comm->nprocs];
  int *displs    = new int[comm->nprocs];
  zlocal         = new double[nlocal];
  zglobal        = new double[natoms];
  deltazs        = new double[natoms];
  
  // Give all processes the number of atoms local to all other processes
  MPI_Allgather(&nlocal,1,MPI_INT,nlocals,1,MPI_INT,world);
  
  // Set up vector of displacements so that when z coordinates are gatheredi,
  //   each processor knows where in the global z vector the z coordinates
  //   from each processor will begin.
  displs[0] = 0;
  for ( i = 1; i < comm->nprocs; i++ ) {
    displs[i] = displs[i-1] + nlocals[i-1];
  }
  
  // Assign z coordinates for local atoms to zlocal, the gather all the zlocals
  //   into a "shared" vector of all atoms' z coordinates, called zglobal.
  for ( i = 0; i < nlocal; i++ ) {
    zlocal[i] = x[i][2];
  }
  MPI_Allgatherv(zlocal,nlocal,MPI_DOUBLE,zglobal,nlocals,displs,MPI_DOUBLE,world);
  
  // The root process sorts the z-positions in order, then uses the differences
  //   between adjacent indices in the vector to obtain a vector of delta-z
  //   values.
  if ( comm->me == 0 ) {
    sortdoubles(zglobal,natoms);
    zl = domain->boxhi[2] - domain->boxlo[2];   // z-length of simulation cell
    i = 0;
    deltazs[i] = zglobal[i] - zglobal[natoms-1];            // Wrap the first
    deltazs[i] = deltazs[i] - zl * floor( deltazs[i] / zl); //  delta-z value.
    sum = deltazs[i];
    sumsq = sum*sum;
    for ( i = 1; i < natoms; i++ ) {
      deltazs[i] = (zglobal[i] - zglobal[i-1]);
      sum += deltazs[i];
      sumsq += deltazs[i]*deltazs[i];
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
  delete[] zlocal;
  delete[] zglobal;
  delete[] deltazs;
  delete[] nlocals;
  delete[] displs;
  
  return scalar;
}


/* ---------------------------------------------------------------------- */
// A function that sorts a vector (z) containing n doubles
void ComputeOPCylDeltaZ::sortdoubles(double *z,bigint n) {
  int i,j;
  double a;
  
  for (i = 0; i < n-1; i++) {
    for (j = i + 1; j < n; j++) {
      if (z[i] > z[j]) {
        a = z[i];
        z[i] = z[j];
        z[j] = a;
      }
    }
  }
}



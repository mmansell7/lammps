/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   This compute was created by J. Matthew ("Matt") Mansell, North
   Caroline State University, August 2019.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_compforce_atom.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCompForceAtom::ComputeCompForceAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute compforce/atom command");

  nmax = 0;

}

/* ---------------------------------------------------------------------- */

ComputeCompForceAtom::~ComputeCompForceAtom()
{
  memory->destroy(compforce);
}


/* ---------------------------------------------------------------------- */

void ComputeCompForceAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCompForceAtom::compute_peratom()
{
  int i;

  // grow local compforce array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(compforce);
    nmax = atom->nmax;
    memory->create(compforce,nmax,"compforce/atom:compforce");
    vector_atom = compforce;
  }

  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info

  int nlocal = atom->nlocal;
  int npair = nlocal;
  if (force->newton) npair += atom->nghost;

  // clear local compforce array

  for (i = 0; i < npair; i++) compforce[i] = 0.0;

  // add in per-atom contributions from each pair force

  if (force->pair) {
    double *cfatom = force->pair->eatom;
    for (i = 0; i < npair; i++) compforce[i] += eatom[i];
  }

  // add in per-atom contributions from relevant fixes
  // always only for owned atoms, not ghost

  if (fixflag && modify->n_thermo_compforce_atom)
    modify->thermo_compforce_atom(nlocal,compforce);

  // communicate ghost compforce between neighbor procs

  if (force->newton || (force->kspace && force->kspace->tip4pflag))
    comm->reverse_comm_compute(this);

  // zero compforce of atoms not in group
  // only do this after comm since ghost contributions must be included

  int *mask = atom->mask;

  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) compforce[i] = 0.0;
}







/* ---------------------------------------------------------------------- */

void ComputeCompForceAtom::pair_contribution()
{
  int i,j,ii,jj,inum,jnum;
  double delx,dely,delz;
  double rsq,eng,cf,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  double cf[4];
  cf[0] = cf[1] = cf[2] = cf[3] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    jlist = firstneigh[i];
    jnum  = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = atom->x[i][0] - atom->x[j][0];
      dely = atom->x[i][1] - atom->x[j][1];
      delz = atom->x[i][2] - atom->x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutsq[atom->type[i]][atom->type[j]]) {
        eng    = pair->single(i,j,atom->type[i],atom->type[j],rsq,factor_coul,factor_lj,fpair);
        
        cf[1] += fabs(delx)*fpair;
        cf[2] += fabs(dely)*fpair;
        cf[3] += fabs(delz)*fpair;
        cf[0] += sqrt(cf[1]*cf[1] + cf[2]*cf[2] + cf[3]*cf[3]);

      }
    }
  }

  double all[4];
  MPI_Allreduce(cf,all,4,MPI_DOUBLE,MPI_SUM,world);
  scalar += all[0];
  vector[0] += all[1];
  vector[1] += all[2];
  vector[2] += all[3];
}

/* ---------------------------------------------------------------------- */







/* ---------------------------------------------------------------------- */

int ComputePEAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = compforce[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePEAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    compforce[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePEAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
















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

#include "compute_compforce_atom.h"
#include <cmath>
#include "atom.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

ComputeCompForceAtom::ComputeCompForceAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), compforce(NULL)
{
  peratom_flag = 1;
  size_peratom_cols = 4;
  comm_reverse = 4;
  timeflag = 1;
  
  if (narg < 3) error->all(FLERR,"Illegal compute compforce/atom command");

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCompForceAtom::~ComputeCompForceAtom()
{
  memory->destroy(compforce);
}

/* ---------------------------------------------------------------------- */

void ComputeCompForceAtom::init()
{
  // need an occasional full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeCompForceAtom::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCompForceAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  
  // grow compforce array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(compforce);
    nmax = atom->nmax;
    memory->create(compforce,nmax,4,"compforce/atom:compforce");
    array_atom = compforce;
  }
  
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,factor_lj,factor_coul,cfx,cfy,cfz,cfmag;
  int *ilist,*jlist,*numneigh,**firstneigh;

  for (i = 0; i < nmax; i++)
    for (j = 0; j < 4; j++)
      compforce[i][j] = 0.0;
      
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < force->pair->cutsq[itype][jtype]) {
        force->pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

        cfx = fabs(delx)*fpair;
        cfy = fabs(dely)*fpair;
        cfz = fabs(delz)*fpair;
        cfmag = sqrt(cfx*cfx + cfy*cfy + cfz*cfz);
        compforce[i][0] += cfmag;
        compforce[i][1] += cfx;
        compforce[i][2] += cfy;
        compforce[i][3] += cfz;

      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCompForceAtom::memory_usage()
{
  double bytes = nmax*4 * sizeof(double);
  return bytes;
}


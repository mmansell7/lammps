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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_pressure_plane.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "bond.h"

using namespace LAMMPS_NS;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// configurational contribution to the normal component of the planar 
// pressure tensor AND tangential configurational component...
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



ComputePressurePlane::ComputePressurePlane(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{

  // Arguments [6]:
  // 0. ID
  // 1. group-ID
  // 2. pressure/plane
  // 3. zmin : minimum z value to include in local pressure analysis
  // 4. zmax : maximum z value to include in local pressure analysis
  // 5. bin_width : width of bins (in r-direction) to include
  //

  if (narg != 6) error->all(FLERR,"Illegal compute pressure/plane command");

  zmin=force->numeric(FLERR,arg[3]);
  zmax=force->numeric(FLERR,arg[4]);
  bin_width=force->numeric(FLERR,arg[5]);

  nbins=(int)((zmax-zmin)/bin_width);

  array_flag=1;
  vector_flag=0;
  extarray=0;
  size_array_cols = 4;  // r, density, PN, PT
  size_array_rows = nbins;

  PNtemp = new double[nbins];
  PNall = new double[nbins];
  PTtemp = new double[nbins];
  PTall = new double[nbins];
  Z = new double[nbins];
  Zkin = new double[nbins];

  density_temp = new double[nbins];
  density_all = new double[nbins];

  memory->create(array,nbins,4,"PN:array");

  nktv2p = force->nktv2p;

  
}

/* ---------------------------------------------------------------------- */

ComputePressurePlane::~ComputePressurePlane()
{
  memory->destroy(array);
  delete [] PNtemp;
  delete [] PNall;
  delete [] PTtemp;
  delete [] PTall;
  delete [] density_temp;
  delete [] density_all;
  delete [] Z;
  delete [] Zkin;
}

/* ---------------------------------------------------------------------- */

void ComputePressurePlane::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"No pair style is defined for compute pressure/plane");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute pressure/plane");

  // precompute some stuff
  for (int iq=0;iq<nbins;iq++)
  {
    Z[iq]=((double)iq+0.5)*bin_width+zmin;
    Zkin[iq]=(((double)iq)+1.0)*bin_width+zmin;
  }

  // for PTAinv, we assume Lx=Ly and only use Lx
  PNAinv=1.0/((domain->boxhi[0]-domain->boxlo[0])*(domain->boxhi[1]-domain->boxlo[1]));
  PTAinv=0.5*PNAinv;
  invVbin=PNAinv/bin_width;

  // need an occasional half neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;

}

/* ---------------------------------------------------------------------- */

void ComputePressurePlane::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
------------------------------------------------------------------------- */

void ComputePressurePlane::compute_array()
{
    invoked_array = update->ntimestep; 

    int ibin;

    // clear pressures
    for (ibin=0;ibin<nbins;ibin++) 
    {
        PNtemp[ibin]=0.0;
        PNall[ibin]=0.0;
        density_temp[ibin]=0.0;
        density_all[ibin]=0.0;
        PTtemp[ibin]=0.0;
        PTall[ibin]=0.0;
    }

    // what processor am I?
    int me;
    MPI_Comm_rank(world,&me);

    int i,j,n,ii,jj,inum,jnum,itype,jtype;
    tagint itag,jtag;
    double xtmp,ytmp,ztmp,delx,dely,delz;
    double rsq,eng,fpair,factor_coul,factor_lj;
    int *ilist,*jlist,*numneigh,**firstneigh;

    double **x = atom->x;
    tagint *tag = atom->tag;
    int *type = atom->type;
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

    // calculate density
    for (i=0;i<nlocal;i++)
    {
        if (x[i][2]>zmax || x[i][2]<zmin) continue; // outside of z limits

        for (j=0;j<nbins;j++) if (x[i][2]<Zkin[j]) break;

        density_temp[j]+=invVbin;
    }
    MPI_Allreduce(density_temp,density_all,nbins,MPI_DOUBLE,MPI_SUM,world);
    for (i=0;i<nbins;i++) array[i][1]=density_all[i];


    // loop over neighbors of my atoms
    // skip if I or J are not in group
    // for newton = 0 and J = ghost atom,
    //   need to insure I,J pair is only output by one proc
    //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()
    // for flag = 0, just count pair interactions within force cutoff
    // for flag = 1, calculate requested output fields

    // PAIR
    double lower_z,upper_z,ft;

    if (force->pair) 
    {
        Pair *pair = force->pair;
        double **cutsq = force->pair->cutsq;

        for (ii = 0; ii < inum; ii++) 
        {
            i = ilist[ii];
            if (!(mask[i] & groupbit)) continue;

            xtmp = x[i][0];
            ytmp = x[i][1];
            ztmp = x[i][2];
            itag = tag[i];
            itype = type[i];
            jlist = firstneigh[i];
            jnum = numneigh[i];

            //r1=x[i][0]*x[i][0]+x[i][1]*x[i][1]+x[i][2]*x[i][2];

            for (jj = 0; jj < jnum; jj++) 
            {
                j = jlist[jj];
                factor_lj = special_lj[sbmask(j)];
                factor_coul = special_coul[sbmask(j)];
                j &= NEIGHMASK;

                if (!(mask[j] & groupbit)) continue;

                // itag = jtag is possible for long cutoffs that include images of self
                // do calculation only on appropriate processor
                if (newton_pair == 0 && j >= nlocal) 
                {
                jtag = tag[j];
                if (itag > jtag) 
                {
                  if ((itag+jtag) % 2 == 0) continue;
                } 
                else if (itag < jtag) 
                {
                  if ((itag+jtag) % 2 == 1) continue;
                }  
                else 
                {
                  if (x[j][2] < ztmp) continue;
                  if (x[j][2] == ztmp) 
                  {
                    if (x[j][1] < ytmp) continue;
                    if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
                  }
                }
                }

                delx = xtmp - x[j][0];
                dely = ytmp - x[j][1];
                delz = ztmp - x[j][2];

                rsq = delx*delx + dely*dely + delz*delz;
                jtype = type[j];
                if (rsq >= cutsq[itype][jtype]) continue;

                eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

                if (x[i][2]<x[j][2])  lower_z=x[i][2];
                else                  lower_z=x[j][2];
                if (lower_z==x[i][2]) upper_z=x[j][2];
                else                  upper_z=x[i][2];

                if (lower_z>zmax || upper_z<zmin) continue; // no contribution

                // pressure contribution
                for (ibin=0;ibin<nbins;ibin++) if (Z[ibin]>lower_z && Z[ibin]<upper_z)
                {
                    PNtemp[ibin]+=fabs(delz)*fpair;
                    //PTtemp[ibin]+=(fabs(delx)+fabs(dely))*fpair;
                    ft=delx*delx+dely*dely;
                    ft/=fabs(delz);
                    ft*=fpair;
                    PTtemp[ibin]+=ft;
                }
            }
        }
    }

    // BONDS 
    if (force->bond)
    {
        int **bondlist = neighbor->bondlist;
        int nbondlist = neighbor->nbondlist;
        int n,type;
        double fbond;

        for (n = 0; n < nbondlist; n++)
        {
            i = bondlist[n][0];
            j = bondlist[n][1];
            type = bondlist[n][2];

            delx = x[i][0] - x[j][0];
            dely = x[i][1] - x[j][1];
            delz = x[i][2] - x[j][2];

            rsq = delx*delx + dely*dely + delz*delz;
            eng = force->bond->single(type,rsq,i,j,fbond);

            if (x[i][2]<x[j][2])  lower_z=x[i][2];
            else                  lower_z=x[j][2];
            if (lower_z==x[i][2]) upper_z=x[j][2];
            else                  upper_z=x[i][2];


            // add forces to pressure contributions
            for (ibin=0;ibin<nbins;ibin++) if (Z[ibin]>lower_z && Z[ibin]<upper_z)
            PNtemp[ibin]+=fabs(delz)*fbond;
            ft=delx*delx+dely*dely;
            ft/=fabs(delz);
            ft*=fbond;
            PTtemp[ibin]+=ft;
        }
    }


    // calculate pressure (force over area)
    for (ibin=0;ibin<nbins;ibin++) 
    {
        PNtemp[ibin]*=PNAinv;
        PTtemp[ibin]*=PTAinv;
    }

    // communicate these values across processors
    MPI_Allreduce(PNtemp,PNall,nbins,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(PTtemp,PTall,nbins,MPI_DOUBLE,MPI_SUM,world);

    // populate array
    for (ibin=0;ibin<nbins;ibin++) 
    {
        array[ibin][0]=Z[ibin];
        array[ibin][2]=PNall[ibin]*nktv2p;
        array[ibin][3]=PTall[ibin]*nktv2p;
    }

}
              

/* ----------------------------------------------------------------------
memory usage of data
------------------------------------------------------------------------- */

double ComputePressurePlane::memory_usage()
{
    double bytes = (1080.0 + 17.0*(double)nbins) * sizeof(double);
    return bytes;
}

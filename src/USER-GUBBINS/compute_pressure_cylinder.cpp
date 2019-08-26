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
#include "compute_pressure_cylinder.h"
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

using namespace LAMMPS_NS;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// configurational contribution to the normal component of the cylindrical
// pressure tensor, PN and tangential configurational component...
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



ComputePressureCylinder::ComputePressureCylinder(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // Arguments [7]:
  // 0. ID
  // 1. group-ID
  // 2. pressure/cylinder
  // 3. zlo : minimum z value to include in local pressure analysis
  // 4. zhi : maximum z value to include in local pressure analysis
  // 5. Rmax : maximum r value to include in analysis
  // 6. bin_width : width of bins (in r-direction) to include
  //

  if (narg != 7) error->all(FLERR,"Illegal compute pressure/cylinder command");

  zlo=force->numeric(FLERR,arg[3]);
  zhi=force->numeric(FLERR,arg[4]);
  Rmax=force->numeric(FLERR,arg[5]);
  bin_width=force->numeric(FLERR,arg[6]);

  nbins=(int)(Rmax/bin_width);

  nzbins=(int)((zhi-zlo)/bin_width);

  array_flag=1;
  vector_flag=0;
  extarray=0;
  size_array_cols = 14; // [0] r
                        // [1] density
                        // [2] Prr_config
                        // [3] Pphiphi_config
                        // [4] Pzz_config, 
                        // [5] Prphi_config
                        // [6] Prz_config
                        // [7] Pphiz_config
                        // [8] Prr_kin
                        // [9] Pphiphi_kin
                        // [10] Pzz_kin
                        // [11] Prphi_kin
                        // [12] Prz_kin
                        // [13] Pphiz_kin
  size_array_rows = nbins;

  Prr_ctemp = new double[nbins];
  Prr_ktemp = new double[nbins];
  Prr_config = new double[nbins];
  Prr_kin = new double[nbins];
  Pphiphi_ctemp = new double[nbins];
  Pphiphi_ktemp = new double[nbins];
  Pphiphi_config = new double[nbins];
  Pphiphi_kin = new double[nbins];
  Pzz_ctemp = new double[nbins];
  Pzz_ktemp = new double[nbins];
  Pzz_config = new double[nbins];
  Pzz_kin = new double[nbins];
  
  Prphi_ctemp = new double[nbins];
  Prphi_ktemp = new double[nbins];
  Prphi_config = new double[nbins];
  Prphi_kin = new double[nbins];
  Prz_ctemp = new double[nbins];
  Prz_ktemp = new double[nbins];
  Prz_config = new double[nbins];
  Prz_kin = new double[nbins];
  Pphiz_ctemp = new double[nbins];
  Pphiz_ktemp = new double[nbins];
  Pphiz_config = new double[nbins];
  Pphiz_kin = new double[nbins];
  
  R  = new double[nbins];
  R2 = new double[nbins];
  PrAinv = new double[nbins];
  PzAinv = new double[nbins];
  Rinv = new double[nbins];
  binz = new double[nzbins];

  R2kin = new double[nbins];
  density_temp = new double[nbins];
  invVbin = new double[nbins];
  density_all = new double[nbins];

  memory->create(array,size_array_rows,size_array_cols,"PN:array");

  nphi=360;
  tangent = new double[nphi];
  ephi_x = new double[nphi];
  ephi_y = new double[nphi];

  nktv2p = force->nktv2p;
  
}

/* ---------------------------------------------------------------------- */

ComputePressureCylinder::~ComputePressureCylinder()
{
  // count all of these for memory usage <------
  memory->destroy(array);
  delete [] R;
  delete [] Rinv;
  delete [] R2;
  delete [] R2kin;
  delete [] invVbin;
  delete [] density_temp;
  delete [] density_all;
  delete [] tangent;
  delete [] ephi_x;
  delete [] ephi_y;
  delete [] Prr_ctemp;
  delete [] Prr_ktemp;
  delete [] Prr_config;
  delete [] Prr_kin;
  delete [] Pzz_ctemp;
  delete [] Pzz_ktemp;
  delete [] Pzz_config;
  delete [] Pzz_kin;
  delete [] Pphiphi_ctemp;
  delete [] Pphiphi_ktemp;
  delete [] Pphiphi_config;
  delete [] Pphiphi_kin;
  delete [] Prphi_ctemp;
  delete [] Prphi_ktemp;
  delete [] Prphi_config;
  delete [] Prphi_kin;
  delete [] Prz_ctemp;
  delete [] Prz_ktemp;
  delete [] Prz_config;
  delete [] Prz_kin;
  delete [] Pphiz_ctemp;
  delete [] Pphiz_ktemp;
  delete [] Pphiz_config;
  delete [] Pphiz_kin;
  delete [] PrAinv;
  delete [] PzAinv;
  delete [] binz;
}

/* ---------------------------------------------------------------------- */

void ComputePressureCylinder::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"No pair style is defined for compute pressure/cylinder");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute pressure/cylinder");

  double phi;

  for (int iphi=0;iphi<nphi;iphi++)
  {
      phi=((double)iphi)*3.14159/180.0;
      tangent[iphi]=tan(phi);
      ephi_x[iphi]=-sin(phi);
      ephi_y[iphi]=cos(phi);
  }

  // precompute some stuff
  for (int iq=0;iq<nbins;iq++)
  {
      R[iq]=((double)iq+0.5)*bin_width;
      Rinv[iq]=1.0/R[iq];
      R2[iq]=R[iq]*R[iq];
      R2kin[iq]=(((double)iq)+1.0)*bin_width;
      R2kin[iq]*=R2kin[iq];
      PrAinv[iq]=1.0/(2.0*3.14159*(zhi-zlo)*R[iq]); // NEW
  }
  PphiAinv=1.0/((zhi-zlo)*bin_width*2.0*(double)nphi); // NEW

  invVbin[0]=1.0/((zhi-zlo)*3.14159*R2kin[0]); // NEW
  PzAinv[0]=1.0/(3.14159*R2kin[0]*((double)nzbins)); // NEW 
  for (int jq=1;jq<nbins;jq++) 
  {
      invVbin[jq]=1.0/((zhi-zlo)*3.14159*(R2kin[jq]-R2kin[jq-1])); // NEW
      PzAinv[jq]=1.0/(3.14159*(R2kin[jq]-R2kin[jq-1])*((double)nzbins)); // NEW
  }

  // need an occasional half neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;

  for (int zzz=0;zzz<nzbins;zzz++) binz[zzz]=(((double)zzz)+0.5)*bin_width+zlo;

}

/* ---------------------------------------------------------------------- */

void ComputePressureCylinder::init_list(int id, NeighList *ptr)
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

void ComputePressureCylinder::compute_array()
{
    invoked_array = update->ntimestep; 

    int ibin;

    // clear pressures
    for (ibin=0;ibin<nbins;ibin++) 
    {
        density_temp[ibin]=0.0; // NEW
        density_all[ibin]=0.0; // NEW
        Prr_ctemp[ibin]=0.0; // NEW
        Prr_ktemp[ibin]=0.0; // NEW
        Prr_config[ibin]=0.0; // NEW
        Prr_kin[ibin]=0.0; // NEW
        Pphiphi_ctemp[ibin]=0.0; // NEW
        Pphiphi_ktemp[ibin]=0.0; // NEW
        Pphiphi_config[ibin]=0.0; // NEW
        Pphiphi_kin[ibin]=0.0; // NEW
        Pzz_ctemp[ibin]=0.0; // NEW
        Pzz_ktemp[ibin]=0.0; // NEW
        Pzz_config[ibin]=0.0; // NEW
        Pzz_kin[ibin]=0.0; // NEW
        Prphi_ctemp[ibin]=0.0; // NEW
        Prphi_ktemp[ibin]=0.0; // NEW
        Prphi_config[ibin]=0.0; // NEW
        Prphi_kin[ibin]=0.0; // NEW
        Prz_ctemp[ibin]=0.0; // NEW
        Prz_ktemp[ibin]=0.0; // NEW
        Prz_config[ibin]=0.0; // NEW
        Prz_kin[ibin]=0.0; // NEW
        Pphiz_ctemp[ibin]=0.0; // NEW
        Pphiz_ktemp[ibin]=0.0; // NEW
        Pphiz_config[ibin]=0.0; // NEW
        Pphiz_kin[ibin]=0.0; // NEW
    }

    // what processor am I?
    int me;
    MPI_Comm_rank(world,&me);

    int i,j,n,ii,jj,inum,jnum,itype,jtype;
    tagint itag,jtag;
    double xtmp,ytmp,ztmp,delx,dely,delz;
    double rsq,eng,fpair,factor_coul,factor_lj;
    int *ilist,*jlist,*numneigh,**firstneigh;
    double temp_phi,vx,vy,vz,vr,vphi;

    double **x = atom->x;
    double *mass = atom->mass;
    double **vtmp = atom->v;
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

    // calculate density and kinetic pressure (by radius)
    double temp_R2;
    for (i=0;i<nlocal;i++) if (x[i][2]<zhi && x[i][2]>zlo)
    {
        temp_R2=x[i][0]*x[i][0]+x[i][1]*x[i][1]; // NEW
        if (temp_R2>R2kin[nbins-1]) continue; // outside of Rmax
        
        temp_phi = atan2(x[i][1],x[i][0]);
        
        vx = vtmp[i][0];
        vy = vtmp[i][1];
        vz = vtmp[i][2];
        
        vr = vx*cos(temp_phi) + vy*sin(temp_phi);
        vphi = -vx*sin(temp_phi) + vy*cos(temp_phi);
        
        for (j=0;j<nbins;j++) if (temp_R2<R2kin[j]) break;
        
        density_temp[j]  += invVbin[j];
        Prr_ktemp[j]     += vr*vr * mass[type[i]] * force->mvv2e * invVbin[j] * nktv2p;
        //printf("Setting Prr_ktemp[%d] += %.5f * %.5g * %.5f * %.5f * %.5f = %.5f\n",
        //          j,vr,vr,mass[type[i]],force->mvv2e,invVbin[j],Prr_ktemp[j]);
        Pphiphi_ktemp[j] += vphi*vphi * mass[type[i]] * force->mvv2e * invVbin[j] * nktv2p;
        Pzz_ktemp[j]     += vz*vz * mass[type[i]] * force->mvv2e * invVbin[j] * nktv2p;
        
        Prphi_ktemp[j]   += vr*vphi * mass[type[i]] * force->mvv2e * invVbin[j] * nktv2p;
        Prz_ktemp[j]     += vr*vz * mass[type[i]] * force->mvv2e * invVbin[j] * nktv2p;
        Pphiz_ktemp[j]   += vphi*vz * mass[type[i]] * force->mvv2e * invVbin[j] * nktv2p;

    }
    MPI_Allreduce(density_temp,density_all,nbins,MPI_DOUBLE,MPI_SUM,world); // NEW
    
    MPI_Allreduce(Prr_ktemp    ,Prr_kin    ,nbins,MPI_DOUBLE,MPI_SUM,world); // NEW
    //printf("Prr_kin mpi_reduced (SUMMED) to %.5f\n",Pr_kin);
    MPI_Allreduce(Pphiphi_ktemp,Pphiphi_kin,nbins,MPI_DOUBLE,MPI_SUM,world); // NEW
    MPI_Allreduce(Pzz_ktemp    ,Pzz_kin    ,nbins,MPI_DOUBLE,MPI_SUM,world); // NEW
    
    MPI_Allreduce(Prphi_ktemp  ,Prphi_kin  ,nbins,MPI_DOUBLE,MPI_SUM,world); // NEW
    MPI_Allreduce(Prz_ktemp    ,Prz_kin    ,nbins,MPI_DOUBLE,MPI_SUM,world); // NEW
    MPI_Allreduce(Pphiz_ktemp  ,Pphiz_kin  ,nbins,MPI_DOUBLE,MPI_SUM,world); // NEW
    
    for (i=0;i<nbins;i++)
    {
        array[i][1]=density_all[i]; // NEW
        
        array[i][8]=Prr_kin[i];
        //printf("Setting array[%d][%d] = %.5f\n",i,5,Pr_kin[i]);
        array[i][9]=Pphiphi_kin[i];
        array[i][10]=Pzz_kin[i];
        
        array[i][11]=Prphi_kin[i];
        array[i][12]=Prz_kin[i];
        array[i][13]=Pphiz_kin[i];
    }

    // loop over neighbors of my atoms
    // skip if I or J are not in group
    // for newton = 0 and J = ghost atom,
    //   need to insure I,J pair is only output by one proc
    //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()
    // for flag = 0, just count pair interactions within force cutoff
    // for flag = 1, calculate requested output fields

    Pair *pair = force->pair;
    double **cutsq = force->pair->cutsq;

    double r1=0.0;
    double r2=0.0;
    double risq,rjsq;
    double ri,rj,rij,fij;
    double A,B,C,Bsq,A2inv,A4,D;
    double alpha1,alpha2,aij;
    double xi,yi,zi,dx,dy,dz;
    double m,xR,yR,zR,fn;
    double alpha,xL,yL,zL,L2,ftphi,ftz;
    double sqrtD;
    double lower_z,upper_z;

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

        r1=x[i][0]*x[i][0]+x[i][1]*x[i][1]; //+x[i][2]*x[i][2];

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


            r2=x[j][0]*x[j][0]+x[j][1]*x[j][1]; //+x[j][2]*x[j][2];

            // ri is smaller of r1 and r2
            if (r2<r1)
            {
                risq=r2;
                rjsq=r1;
                xi=x[j][0];
                yi=x[j][1];
                zi=x[j][2];
                dx=x[i][0]-x[j][0];
                dy=x[i][1]-x[j][1];
                dz=x[i][2]-x[j][2];
            }
            else
            {
                risq=r1;
                rjsq=r2;
                xi=x[i][0];
                yi=x[i][1];
                zi=x[i][2];
                dx=x[j][0]-x[i][0];
                dy=x[j][1]-x[i][1];
                dz=x[j][2]-x[i][2];
            }

            rsq = delx*delx + dely*dely + delz*delz;
            jtype = type[j];
            if (rsq >= cutsq[itype][jtype]) continue;

            eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

            A=dx*dx+dy*dy; //+dz*dz;
            B=2.0*(xi*dx+yi*dy); //+zi*dz);

            // normal pressure contribution PN ====== PN ====== PN ====== PN ====== PN ====== PN ====== PN ====== PN
            for (ibin=0;ibin<nbins;ibin++)
            {
                // completely inside of R
                if (rjsq<R2[ibin]) continue;

                C=risq-R2[ibin];
                D=B*B-4.0*A*C;

                // completely oustide R
                if (D<0.0) continue;

                sqrtD=sqrt(D);
                alpha1=0.5*(-B+sqrtD)/A;
                alpha2=0.5*(-B-sqrtD)/A;

                if (alpha1>0.0 && alpha1<1.0)
                {
                    zR=zi+alpha1*dz;
                    if (zR<zhi && zR>zlo)
                    {
                        xR=xi+alpha1*dx;
                        yR=yi+alpha1*dy;
                        // fpair is force / separation
                        // So fn is force * distance
                        fn=fpair*fabs(xR*dx+yR*dy);

                        Prr_ctemp[ibin]+=fn;

                        //ftz=fabs(dz)*fpair;
                    }
                }
                if (alpha2>0.0 && alpha2<1.0)
                {
                    zR=zi+alpha2*dz;
                    if (zR<zhi && zR>zlo)
                    {
                        xR=xi+alpha2*dx;
                        yR=yi+alpha2*dy;
                        fn=fpair*fabs(xR*dx+yR*dy);

                        Prr_ctemp[ibin]+=fn;

                        //ftz=fabs(dz)*fpair;
                    }
                }
            }

            // tangential pressure contribution PT ====== PT ====== PT ====== PT ====== PT ====== PT ====== PT ====== PT
            for (int iphi=0;iphi<nphi;iphi++)
            {
                // in a VERY rare case, this could be singular, but C is robust
                alpha=(yi-xi*tangent[iphi])/(dx*tangent[iphi]-dy);

                // no intersection with phi surface
                if (alpha>=1.0 || alpha<=0.0) continue;

                // no contribution (outside of averaging region)
                zL=zi+alpha*dz;
                if (zL>zhi || zL<zlo) continue;

                xL=xi+alpha*dx;
                yL=yi+alpha*dy;

                L2=xL*xL+yL*yL; // NEW

                // no intersection (outside of Rmax)
                if (L2>R2kin[nbins-1]) continue;

                ftphi=fabs(dx*ephi_x[iphi]+dy*ephi_y[iphi])*fpair; // NEW

                // add to appropriate bin
                for (ibin=0;ibin<nbins;ibin++) if (L2<R2kin[ibin])
                {
                    Pphiphi_ctemp[ibin]+=ftphi; // NEW
                    break;
                }
            }

            // z pressure contribution
            // use discrete bins instead of a smarter way, for now
            for (int zbin=0;zbin<nzbins;zbin++)
            {
                // check if interaction contributes
                if (x[i][2]>binz[zbin] && x[j][2]>binz[zbin]) continue;
                if (x[i][2]<binz[zbin] && x[j][2]<binz[zbin]) continue;

                alpha=(binz[zbin]-zi)/dz;

                xL=xi+alpha*dx;
                yL=yi+alpha*dy;

                L2=xL*xL+yL*yL; 

                if (L2>R2kin[nbins-1]) continue;

                ftz=fabs(dz)*fpair;

                // add to appropriate bin
                for (ibin=0;ibin<nbins;ibin++) if (L2<R2kin[ibin])
                {
                    Pzz_ctemp[ibin]+=ftz; // NEW
                    break;
                }
            }
        }
    }

    // calculate pressure (force over area)
    for (ibin=0;ibin<nbins;ibin++) 
    {
        // force*distance * 1/area * 1/length = force / area
        Prr_ctemp[ibin]*=PrAinv[ibin]*Rinv[ibin];
        // force * 1/area = force / area
        Pphiphi_ctemp[ibin]*=PphiAinv;
        // force * 1/area = force / area
        Pzz_ctemp[ibin]*=PzAinv[ibin];
    }

    // communicate these values across processors
    MPI_Allreduce(Prr_ctemp,Prr_config,nbins,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(Pphiphi_ctemp,Pphiphi_config,nbins,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(Pzz_ctemp,Pzz_config,nbins,MPI_DOUBLE,MPI_SUM,world);

    MPI_Allreduce(Prphi_ctemp,Prphi_config,nbins,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(Prz_ctemp  ,Prz_config  ,nbins,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(Pphiz_ctemp,Pphiz_config,nbins,MPI_DOUBLE,MPI_SUM,world);

    // populate array
    for (ibin=0;ibin<nbins;ibin++) 
    {
        array[ibin][0]=R[ibin];
      // nktv2p converts nkt/v to pressure units
      // This is required because all calculations above were done in 
      // units of n, k, T, and v (l^3), but the output should be in 
      // pressure units.  In LAMMPS, units of nkt/v are not necessarily
      // equivalent to pressure units.  See LAMMPS "units" command.
        array[ibin][2]=Prr_config[ibin]*nktv2p;
        array[ibin][3]=Pphiphi_config[ibin]*nktv2p;
        array[ibin][4]=Pzz_config[ibin]*nktv2p;
        array[ibin][5]=Prphi_config[ibin]*nktv2p;
        array[ibin][6]=Prz_config[ibin]*nktv2p;
        array[ibin][7]=Pphiz_config[ibin]*nktv2p;
    }

}
              

/* ----------------------------------------------------------------------
memory usage of data
------------------------------------------------------------------------- */

double ComputePressureCylinder::memory_usage()
{
    double bytes = (3.0*(double)nphi + 16.0*(double)nbins+5.0*(double)nbins) * sizeof(double);
    return bytes;
}


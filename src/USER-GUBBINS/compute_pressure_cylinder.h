/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(pressure/cylinder,ComputePressureCylinder)

#else

#ifndef LMP_COMPUTE_PRESSURE_CYLINDER
#define LMP_COMPUTE_PRESSURE_CYLINDER

#include "compute.h"

namespace LAMMPS_NS {

class ComputePressureCylinder : public Compute {
 public:
  ComputePressureCylinder(class LAMMPS *, int, char **);
  ~ComputePressureCylinder();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  double memory_usage();

 private:
  int nbins,nphi,nzbins;
  double *Prr_ctemp    ,*Prr_ktemp    ,*Prr_config    ,*Prr_kin    ;
  double *Pphiphi_ctemp,*Pphiphi_ktemp,*Pphiphi_config,*Pphiphi_kin;
  double *Pzz_ctemp    ,*Pzz_ktemp    ,*Pzz_config    ,*Pzz_kin    ;
  double *Prphi_ctemp,*Prphi_ktemp,*Prphi_config,*Prphi_kin;
  double *Prz_ctemp  ,*Prz_ktemp  ,*Prz_config  ,*Prz_kin  ;
  double *Pphiz_ctemp,*Pphiz_ktemp,*Pphiz_config,*Pphiz_kin;
  double *R,*Rinv,*R2,*PrAinv,*PzAinv,PphiAinv;
  double Rmax,bin_width,nktv2p;
  double *R2kin,*density_temp,*invVbin,*density_all;
  double *tangent,*ephi_x,*ephi_y;
  double *binz;

  double zlo,zhi;

  class NeighList *list;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No pair style is defined for compute pressure/cylinder

Self-explanatory.

E: Pair style does not support compute pressure/cylinder

The pair style does not have a single() function, so it can
not be invoked by compute pressure/cylinder.

*/


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

ComputeStyle(pressure/plane,ComputePressurePlane)

#else

#ifndef LMP_COMPUTE_PRESSURE_PLANE
#define LMP_COMPUTE_PRESSURE_PLANE

#include "compute.h"

namespace LAMMPS_NS {

class ComputePressurePlane : public Compute {
 public:
  ComputePressurePlane(class LAMMPS *, int, char **);
  ~ComputePressurePlane();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  double memory_usage();

 private:
  int nbins,nphi;
  double *PNtemp,*PNall,*PTtemp,*PTall;
  double zmin,zmax,bin_width,nktv2p,PTAinv,PNAinv,invVbin;
  double *density_all,*density_temp,*Zkin,*Z;

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

E: No pair style is defined for compute pressure/sphere

Self-explanatory.

E: Pair style does not support compute pressure/sphere

The pair style does not have a single() function, so it can
not be invoked by compute pressure/sphere.

*/


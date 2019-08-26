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

ComputeStyle(opcyl/deltazt,ComputeOPCylDeltaZT)

#else

#ifndef LMP_COMPUTE_OPCYL_DELTAZT
#define LMP_COMPUTE_OPCYL_DELTAZT

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOPCylDeltaZT : public Compute {
 public:
  ComputeOPCylDeltaZT(class LAMMPS *, int, char **);
  virtual ~ComputeOPCylDeltaZT();
  void init();
  void init_list(int, class NeighList *);
  virtual double compute_scalar();

 private:
  void sortdoubles(double *, bigint);
  double *ztlocal,*ztglobal;
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

W: More than one compute opcyl/deltazt

It is not efficient to use compute deltazt/atom more than once.

*/

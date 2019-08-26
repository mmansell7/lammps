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

ComputeStyle(opcyl/deltaz,ComputeOPCylDeltaZ)

#else

#ifndef LMP_COMPUTE_OPCYL_DELTAZ
#define LMP_COMPUTE_OPCYL_DELTAZ

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOPCylDeltaZ : public Compute {
 public:
  ComputeOPCylDeltaZ(class LAMMPS *, int, char **);
  virtual ~ComputeOPCylDeltaZ();
  void init();
  void init_list(int, class NeighList *);
  virtual double compute_scalar();

 private:
  void sortdoubles(double *, bigint);
  double *zlocal,*zglobal;
  double *deltazs;
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

W: More than one compute opcyl/deltaz

It is not efficient to use compute deltaz/atom more than once.

*/

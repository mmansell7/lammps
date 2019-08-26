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

ComputeStyle(opcyl/deltatz,ComputeOPCylDeltaTZ)

#else

#ifndef LMP_COMPUTE_OPCYL_DELTATZ
#define LMP_COMPUTE_OPCYL_DELTATZ

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOPCylDeltaTZ : public Compute {
 public:
  ComputeOPCylDeltaTZ(class LAMMPS *, int, char **);
  virtual ~ComputeOPCylDeltaTZ();
  void init();
  void init_list(int, class NeighList *);
  virtual double compute_scalar();

 private:
  void sortdoubles(double *, bigint);
  double *tzlocal,*tzglobal;
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

W: More than one compute opcyl/deltatz

It is not efficient to use compute deltatz/atom more than once.

*/

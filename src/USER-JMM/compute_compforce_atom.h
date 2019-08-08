/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file created by J. Matthew ("Matt") Mansell, North Carolina State
   Univerity, August 2019.

------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(compforce/atom,ComputeCompForceAtom)

#else

#ifndef LMP_COMPUTE_COMPFORCE_ATOM_H
#define LMP_COMPUTE_COMPFORCE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCompForceAtom : public Compute {
 public:
  ComputeCompForceAtom(class LAMMPS *, int, char **);
  ~ComputeCompForceAtom();
  void init() {}
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int kspaceflag,fixflag;
  int nmax;
  double *energy;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Per-atom compressive force was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/

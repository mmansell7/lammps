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

ComputeStyle(deltathetaz/atom,ComputeDeltaThetaZAtom)

#else

#ifndef LMP_COMPUTE_DELTATHETAZ_ATOM_H
#define LMP_COMPUTE_DELTATHETAZ_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDeltaThetaZAtom : public Compute {
 public:
  ComputeDeltaThetaZAtom(class LAMMPS *, int, char **);
  ~ComputeDeltaThetaZAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();
  enum {NONE,CUTOFF,ORIENT};

 private:
  int nmax,ncol;
  double cutsq;
  class NeighList *list;

  int *typelo,*typehi;
  double *cvec,*zpmins;
  double **carray;

  // class ComputeOrientOrderAtom *c_orientorder;
  // char *id_orientorder;
  double threshold;
  double **normv;
  int cstyle,nqlist,l;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute deltathetaz/atom requires a pair style be defined

Self-explanatory.

E: Compute deltathetaz/atom cutoff is longer than pairwise cutoff

Cannot identify neighbors at distances longer than the pair cutoff,
since those atoms are not in the neighbor list.

W: More than one compute deltathetaz/atom

It is not efficient to use compute deltathetaz/atom more than once.

*/

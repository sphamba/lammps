// Copyright Son Pham-Ba 2021

#ifdef PAIR_CLASS

PairStyle(gran/cohesive/history,PairGranCohesiveHistory)

#else

#ifndef LMP_PAIR_GRAN_COHESIVE_HISTORY_H
#define LMP_PAIR_GRAN_COHESIVE_HISTORY_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

  class PairGranCohesiveHistory : public PairGranHookeHistory {
  public:
    PairGranCohesiveHistory(class LAMMPS *);
    ~PairGranCohesiveHistory();
    void compute(int, int);
    void settings(int, char **);
    void coeff(int, char **);
    double single(int, int, int, int, double, double, double, double &);
    // double init_one(int, int);
  protected:
    void allocate();

  private: // per-type coefficients, set in pair coeff command
    // pairstyle parameters
    double **E;
    double **nu;
    double **sigma_n_max;
    double **sigma_t_max;
    double **gamma;
    double **scaling;

    // useful quantities computed directly from params
    double **eps_e;
    double **delta_f;
    double **d_c;
    double **c_factor;
  };

}

#endif // PAIR_CLASS
#endif // LMP_PAIR_GRAN_COHESIVE_HISTORY_H

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/

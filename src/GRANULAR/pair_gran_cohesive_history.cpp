// Copyright Son Phamb-Ba 2021

#include "pair_gran_cohesive_history.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "neighbor.h"
#include "modify.h"
#include "neigh_list.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranCohesiveHistory::PairGranCohesiveHistory(LAMMPS *lmp) :
  PairGranHookeHistory(lmp) {}

/* ---------------------------------------------------------------------- */

PairGranCohesiveHistory::~PairGranCohesiveHistory()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(E);
    memory->destroy(nu);
    memory->destroy(sigma_n_max);
    memory->destroy(sigma_t_max);
    memory->destroy(gamma);
    memory->destroy(scaling);
    memory->destroy(eps_e);
    memory->destroy(delta_f);
    memory->destroy(d_c);
    memory->destroy(c_factor);
  }

  // destructor of PairGranHookeHistory is called after
}

/* ---------------------------------------------------------------------- */

void PairGranCohesiveHistory::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  ev_init(eflag,vflag);

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firstshear = fix_history->firstvalue;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      jtype = type[j];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;
      double delta_e = radsum * eps_e[itype][jtype]; // elastic limit
      double delta_f_adj = delta_e * pow(radsum / d_c[itype][jtype], scaling[itype][jtype]); // adjusted for large radii
      delta_f_adj = MAX(delta_f[itype][jtype], delta_f_adj);
      double radsum_adh = radsum + delta_f_adj;

      if (rsq >= radsum_adh*radsum_adh) {

        // unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[3*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        wr1 = (radi*omega[i][0] + radj*omega[j][0]) * rinv;
        wr2 = (radi*omega[i][1] + radj*omega[j][1]) * rinv;
        wr3 = (radi*omega[i][2] + radj*omega[j][2]) * rinv;

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        mi = rmass[i];
        mj = rmass[j];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
          if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
        }

        meff = mi*mj / (mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

	// Custom force here!
	// radsum is the equilibrium distance
	// double R_eq = sqrt(radi * radj); // equivalent single particle radius
	double R_eq = MIN(radi, radj); // equivalent single particle radius
	double A0 = 2.0 * sqrt(3.0) * R_eq * R_eq; // equivalent contact area, assuming hexagonal packing
	double An = A0 / sqrt(6.0) / (1.0 - 2.0 * nu[itype][jtype]); // to balance between k_e and k_t, to have correct E and nu
	double At = An * (1.0 - 4.0 * nu[itype][jtype]) / (1.0 + nu[itype][jtype]);

	double delta = radsum - r; // normal interpenetration
	double rescale_delta = MIN(1.0, 1.0 + delta / delta_f_adj);
	// to rescale (decrease) forces when 0 < delta < delta_f_adj

	double k_ = E[itype][jtype] / radsum;
	double k_e = An * k_; // stiffness elastic domain
	double k_t = At * k_; // tangential stiffness elastic domain
	double k_f = An * sigma_n_max[itype][jtype] / (delta_f_adj - delta_e); // stiffness fracture domain

	double c_n = c_factor[itype][jtype] * sqrt(meff * k_e); // damping normal direction
	double c_t = c_factor[itype][jtype] * sqrt(meff * k_t); // damping tangential direction

	double F_normal = delta >= -delta_e ? k_e * delta : -k_f * (delta + delta_f_adj);
	// delta_f is adjusted to take care of the case when delta_f < delta_e and conserve sigma_n_max

        // normal forces = Hookian contact + normal velocity damping

        damp = c_n*vnnr*rsqinv;
	ccel = F_normal * rinv - damp; // Normal force in w direction = ccel * delw

        // relative velocities (taking into account rotations and distance between particles?)

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;
        shear = &allshear[3*jj];

        if (shearupdate) {
          shear[0] += vtr1*dt;
          shear[1] += vtr2*dt;
          shear[2] += vtr3*dt;
        }
        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
                      shear[2]*shear[2]);

        // rotate shear displacements (remove normal component)

        rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
        rsht *= rsqinv;
        if (shearupdate) {
          shear[0] -= rsht*delx;
          shear[1] -= rsht*dely;
          shear[2] -= rsht*delz;
        }

        // tangential forces = shear + tangential velocity damping

        fs1 = - (k_t*shear[0] + c_t*vtr1);
        fs2 = - (k_t*shear[1] + c_t*vtr2);
        fs3 = - (k_t*shear[2] + c_t*vtr3);

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = At * sigma_t_max[itype][jtype] * rescale_delta; // rescale sigma_t_max according to delta

        if (fs > fn) {
          if (shrmag != 0.0) {
            shear[0] = (fn/fs) * (shear[0] + c_t*vtr1/k_t) - // I don't know why those lines are here,
              c_t*vtr1/k_t;                                  // and not in "single()", but they are
            shear[1] = (fn/fs) * (shear[1] + c_t*vtr2/k_t) - // present in pair_gran_hooke_history.
              c_t*vtr2/k_t;
            shear[2] = (fn/fs) * (shear[2] + c_t*vtr3/k_t) -
              c_t*vtr3/k_t;
            fs1 *= fn/fs;
            fs2 *= fn/fs;
            fs3 *= fn/fs;
          } else fs1 = fs2 = fs3 = 0.0;
        }

        // forces & torques

	// ccel: normal force magnitude (why this name???)
        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);
        torque[i][0] -= radi*tor1;
        torque[i][1] -= radi*tor2;
        torque[i][2] -= radi*tor3;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= radj*tor1;
          torque[j][1] -= radj*tor2;
          torque[j][2] -= radj*tor3;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairGranCohesiveHistory::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(E,           n+1, n+1, "pair:E");
  memory->create(nu,          n+1, n+1, "pair:nu");
  memory->create(sigma_n_max, n+1, n+1, "pair:sigma_n_max");
  memory->create(sigma_t_max, n+1, n+1, "pair:sigma_t_max");
  memory->create(gamma,       n+1, n+1, "pair:gamma");
  memory->create(scaling,     n+1, n+1, "pair:scaling");
  memory->create(eps_e,       n+1, n+1, "pair:eps_e");
  memory->create(delta_f,     n+1, n+1, "pair:delta_f");
  memory->create(d_c,         n+1, n+1, "pair:d_c");
  memory->create(c_factor,    n+1, n+1, "pair:c_factor");

  onerad_dynamic = new double[n+1];
  onerad_frozen  = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen  = new double[n+1];
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairGranCohesiveHistory::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairGranCohesiveHistory::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      E[i][j]           = utils::numeric(FLERR,arg[2],false,lmp);
      nu[i][j]          = utils::numeric(FLERR,arg[3],false,lmp);
      sigma_n_max[i][j] = utils::numeric(FLERR,arg[4],false,lmp);
      sigma_t_max[i][j] = utils::numeric(FLERR,arg[5],false,lmp);
      gamma[i][j]       = utils::numeric(FLERR,arg[6],false,lmp);
      scaling[i][j]     =-utils::numeric(FLERR,arg[7],false,lmp);
      double restit     = utils::numeric(FLERR,arg[8],false,lmp);

      if (E[i][j] < 0.0 || sigma_n_max[i][j] < 0.0 || sigma_t_max[i][j] < 0.0 || gamma[i][j] < 0.0 ||
	  scaling[i][j] > 0.0)
	error->all(FLERR,"Incorrect args for pair coefficients");

      if (nu[i][j] <= -1.0 || nu[i][j] > 0.25)
	error->all(FLERR,"Incorrect args for pair coefficients: nu must be between -1 and 0.25");

      if (scaling[i][j] < -1.0)
	error->all(FLERR,"Incorrect args for pair coefficients: scaling > 1 not implemented");

      // convert from pressure units to force/distance^2

      E[i][j] /= force->nktv2p;
      sigma_n_max[i][j] /= force->nktv2p;
      sigma_t_max[i][j] /= force->nktv2p;
        
      // compute useful quantities
      eps_e[i][j] = sigma_n_max[i][j] / E[i][j];
      delta_f[i][j] = 4.0 * gamma[i][j] / sigma_n_max[i][j]; // fracture delta
      d_c[i][j] = 4.0 * gamma[i][j] * E[i][j] / (sigma_n_max[i][j] * sigma_n_max[i][j]); // critical diameter
      c_factor[i][j] = 2.0 * (1.0 - restit) / M_PI;

      E[j][i] = E[i][j];
      nu[j][i] = nu[i][j];
      sigma_n_max[j][i] = sigma_n_max[i][j];
      sigma_t_max[j][i] = sigma_t_max[i][j];
      gamma[j][i] = gamma[i][j];
      scaling[j][i] = scaling[i][j];
      eps_e[j][i] = eps_e[i][j];
      delta_f[j][i] = delta_f[i][j];
      d_c[j][i] = d_c[i][j];
      c_factor[j][i] = c_factor[i][j];

      setflag[i][j] = 1;
      count++;
    }
  }

  xmu = 0.0;
  dampflag = 1;

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */
/*
double PairGranCohesiveHistory::init_one(int i, int j)
{
  if (!allocated) allocate();

  // cutoff = sum of max I,J radii for
  // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

  double cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
  cutoff = MAX(cutoff,maxrad_frozen[i]+maxrad_dynamic[j]);
  cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]);
  return cutoff + 1.0e-9;
}
*/
/* ---------------------------------------------------------------------- */

double PairGranCohesiveHistory::single(int i, int j, int itype, int jtype,
				       double rsq,
				       double /*factor_coul*/, double /*factor_lj*/,
				       double &fforce)
{
  double radi,radj,radsum;
  double r,rinv,rsqinv,delx,dely,delz;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double mi,mj,meff,damp,ccel;
  double vtr1,vtr2,vtr3,vrel,shrmag,rsht;
  double fs1,fs2,fs3,fs,fn;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;
  double delta_e = radsum * eps_e[itype][jtype]; // elastic limit
  double delta_f_adj = delta_e * pow(radsum / d_c[itype][jtype], scaling[itype][jtype]); // adjusted for large radii
  delta_f_adj = MAX(delta_f[itype][jtype], delta_f_adj);
  double radsum_adh = radsum + delta_f_adj;

  if (rsq >= radsum_adh*radsum_adh) { // rsq (argument) radius squared
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  double **v = atom->v;
  vr1 = v[i][0] - v[j][0];
  vr2 = v[i][1] - v[j][1];
  vr3 = v[i][2] - v[j][2];

  // normal component

  double **x = atom->x;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  double **omega = atom->omega;
  wr1 = (radi*omega[i][0] + radj*omega[j][0]) * rinv;
  wr2 = (radi*omega[i][1] + radj*omega[j][1]) * rinv;
  wr3 = (radi*omega[i][2] + radj*omega[j][2]) * rinv;

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle

  double *rmass = atom->rmass;
  int *mask = atom->mask;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    // NOTE: insure mass_rigid is current for owned+ghost atoms?
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi*mj / (mi+mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  // Custom force here!
  // radsum is the equilibrium distance
  // double R_eq = sqrt(radi * radj); // equivalent single particle radius
  double R_eq = MIN(radi, radj); // equivalent single particle radius
  double A0 = 2.0 * sqrt(3.0) * R_eq * R_eq; // equivalent contact area, assuming hexagonal packing
  double An = A0 / sqrt(6.0) / (1.0 - 2.0 * nu[itype][jtype]); // to balance between k_e and k_t, to have correct E and nu
  double At = An * (1.0 - 4.0 * nu[itype][jtype]) / (1.0 + nu[itype][jtype]);

  double delta = radsum - r; // normal interpenetration
  double rescale_delta = MIN(1.0, 1.0 + delta / delta_f_adj);
  // to rescale (decrease) forces when 0 < delta < delta_f_adj

  double k_ = E[itype][jtype] / radsum;
  double k_e = An * k_; // stiffness elastic domain
  double k_t = At * k_; // tangential stiffness elastic domain
  double k_f = An * sigma_n_max[itype][jtype] / (delta_f_adj - delta_e); // stiffness fracture domain

  double c_n = c_factor[itype][jtype] * sqrt(meff * k_e); // damping normal direction
  double c_t = c_factor[itype][jtype] * sqrt(meff * k_t); // damping tangential direction

  double F_normal = delta >= -delta_e ? k_e * delta : -k_f * (delta + delta_f_adj);
  // delta_f is adjusted to take care of the case when delta_f < delta_e and conserve sigma_n_max

  // normal forces = Hookian contact + normal velocity damping

  damp = c_n*vnnr*rsqinv;
  ccel = F_normal * rinv - damp; // Normal force in w direction = ccel * delw

  // relative velocities

  vtr1 = vt1 - (delz*wr2-dely*wr3);
  vtr2 = vt2 - (delx*wr3-delz*wr1);
  vtr3 = vt3 - (dely*wr1-delx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects
  // neighprev = index of found neigh on previous call
  // search entire jnum list of neighbors of I for neighbor J
  // start from neighprev, since will typically be next neighbor
  // reset neighprev to 0 as necessary

  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  double *allshear = fix_history->firstvalue[i];

  for (int jj = 0; jj < jnum; jj++) {
    neighprev++;
    if (neighprev >= jnum) neighprev = 0;
    if (jlist[neighprev] == j) break;
  }

  double *shear = &allshear[3*neighprev];
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
                shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
  rsht *= rsqinv;

  // tangential forces = shear + tangential velocity damping

  fs1 = - (kt*shear[0] + c_t*vtr1);
  fs2 = - (kt*shear[1] + c_t*vtr2);
  fs3 = - (kt*shear[2] + c_t*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = At * sigma_t_max[itype][jtype] * rescale_delta; // rescale sigma_t_max according to delta

  if (fs > fn) {
    if (shrmag != 0.0) {
      fs1 *= fn/fs;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
      fs *= fn/fs;
    } else fs1 = fs2 = fs3 = fs = 0.0;
  }

  // set force and return no energy

  fforce = ccel;

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vt1;
  svector[8] = vt2;
  svector[9] = vt3;

  return 0.0;
}

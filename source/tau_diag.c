/** ************************************************************************* */
/**
 * @file     tau_diag.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    The main resting place for all functions related to providing
 *           optical depth diagnostic.
 *
 * ************************************************************************** */

#include <stdio.h>

#include "atomic.h"
#include "python.h"

/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]
 *
 * @return
 *
 * @details
 *
 * ************************************************************************** */

void
init_tau_diag_photon (PhotPtr pout)
{
  pout->x[0] = geo.rstar;
  pout->x[1] = 0;
  pout->x[2] = 0;

  pout->lmn[0] = 0;
  pout->lmn[1] = 0;
  pout->lmn[2] = 0;

  pout->w = geo.f_tot;
  pout->w_orig = geo.f_tot;
  pout->freq = 3.2e15;
  pout->freq_orig = 3.2e15;
  pout->tau = 0;

  pout->istat = P_INWIND;
  pout->nres = -1;  // electron scattering
  pout->nrscat = 0;
  pout->nscat = 0;
  pout->nnscat = 0;

  pout->grid = 0;
  pout->origin = PTYPE_DISK;
  pout->origin_orig = PTYPE_DISK;

  pout->nnscat = 0;
  pout->np = 0;

  pout->path = 0;
}

/* ************************************************************************* */
/**
 * @brief     The main controlling function for the optical depth diagnostics.
 *
 * @return    void
 *
 * @details
 *
 * This is the main function which will control the procedure for calculating
 * various diagnostic numbers for the optical depth's experienced in the current
 * model. Namely, this function aims to show the total integrated optical depth
 * to each observer angle using the following optical depths:
 *
 *  - Rosseland mean
 *  - Planck mean
 *  - Lymann edge
 *  - Balmer edge
 *  - Average freq for the photon bands
 *
 * The aim of these diagnostic numbers it to provide some sort of metric on the
 * optical thickness of the current model.
 *
 * ************************************************************************** */

void
tau_diag (void)
{
  int ispec;
  struct photon p, pp;
  double inclination;
  WindPtr w;

  w = wmain;

  Log ("\nThe optical depth along the observer LOS's:\n");

  init_tau_diag_photon (&p);

  // MSPEC here is used as this is where the actual observer angles start
  for (ispec = 0; ispec < nspectra - MSPEC; ispec++)
  {
    stuff_phot (&p, &pp);
    stuff_v (xxspec[ispec].lmn, pp.lmn);

    inclination = geo.angle[ispec];
    Log ("geo.angle[%i] = %f\n", ispec, inclination);

    extract_one (w, &p, PTYPE_DISK, ispec + MSPEC);
    Log ("p->tau = %e\n", &p.tau);
  }

}

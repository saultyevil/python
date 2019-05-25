/** ************************************************************************* */
/**
 * @file     mean_opacities.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    Contains functions used to calculate mean opaicities, such as the
 *           Rosseland mean.
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"


#define DEFAULT_DOMAIN 0

double *wind_rosseland_opac;
double *wind_planck_opac;

/* ************************************************************************* */
/**
 * @brief           Allocate memory for the global opacity arrays.
 *
 * @return          int               EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * ************************************************************************** */

int
allocate_opacity_arrays (int ncells)
{
  wind_rosseland_opac = calloc (ncells, sizeof (*wind_rosseland_opac));
  if (wind_rosseland_opac == NULL)
  {
    Error ("%s:%i:%s: cannot allocate %d bytes for wind_rosseland_opac\n", __FILE__, __func__, __LINE__, ncells *
           sizeof (*wind_rosseland_opac));
    return EXIT_FAILURE;
  }

  Log_silent ("Allocated %d bytes for wind_rosseland_opac\n", ncells * sizeof (*wind_rosseland_opac));

  wind_planck_opac = calloc (ncells, sizeof (*wind_planck_opac));
  if (wind_planck_opac == NULL)
  {
    Error ("%s:%i:%s: cannot allocate %d bytes for wind_planck_opac\n", __FILE__, __func__, __LINE__, ncells *
           sizeof (*wind_planck_opac));
    return EXIT_FAILURE;
  }

  Log_silent ("Allocated %d bytes for wind_planck_opac\n", ncells * sizeof (*wind_planck_opac));

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief     Free memory for the global opacity arrays.
 *
 * @return    void
 *
 * @details
 *
 * ************************************************************************** */

void
free_opacity_arrays (void)
{
  free (wind_rosseland_opac);
  free (wind_planck_opac);
}


/* ************************************************************************* */
/**
 * @brief     Computes the temperature derivative of the Planck function.
 *
 * @param[in] double nu             The frequency of the black body emission
 * @param[in] double temperature    The temperature of the black body
 *
 * @returns   double dB_dT          The temperature derivative of the Planck
 *                                  function wrt to temperature for a given
 *                                  frequency and temperature
 *
 * @details
 *
 * I take some shortcuts with algebra for dB/dT hence the final expression is a
 * bit odd looking.
 *
 * ************************************************************************** */

double ross_int_temp;

double
dB_nu_dT (double nu)
{
  double x;
  double bnu;
  double dB_dT;

  x = H_OVER_K * nu / ross_int_temp;
  bnu = 2 * H * pow (nu, 3) / pow (C, 2) * 1 / (exp (x) - 1);
  dB_dT = bnu * x / (ross_int_temp * (1 - exp (-x)));

  return dB_dT;
}

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
set_ross_integ_temp (double temperature)
{
  ross_int_temp = temperature;
}

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
compute_rosseland_opacity (void)
{
}


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
compute_planck_opacity (void)
{
}

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

int
create_test_photon (PhotPtr pout, double nu, double x[3])
{
  if (nu < 0)
  {
    Error ("%s:%s:%i: photon can't be created with negative nu\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  pout->x[0] = x[0];
  pout->x[1] = x[1];
  pout->x[2] = x[2];
  pout->lmn[0] = 1.0;
  pout->lmn[1] = 0.0;
  pout->lmn[2] = 0.0;
  pout->freq = pout->freq_orig = nu;
  pout->tau = 0.0;
  pout->istat = P_INWIND;
  pout->origin = pout->origin_orig = PTYPE_WIND;

  return EXIT_SUCCESS;
}

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

int
init_mean_opacities (void)
{
  int i;
  int nx;
  int nz;
  int ncells;
  int wind_start_index;
  int wind_stop_index;

  double ds;
  double *cell_x;
  double *next_cell_x;
  double opac_freq;

  struct photon popac;

  Log_silent ("Initialising mean opacities for all wind cells\n");

  /*
   * We need to know how many cells we are going to be calculating the mean
   * opacities for...
   */

  nx = zdom[DEFAULT_DOMAIN].ndim;
  nz = zdom[DEFAULT_DOMAIN].mdim;
  wind_start_index = zdom[DEFAULT_DOMAIN].nstart;
  wind_stop_index = zdom[DEFAULT_DOMAIN].nstop;
  ncells = nx * nz;


  if (allocate_opacity_arrays (ncells))
    return EXIT_FAILURE;

  /*
   * Loop over all cells in the simulation... and fire multiple test photons to
   * calculate the frequency dependent opacity for a range of frequency photons
   * TODO: check that i should actually reach nstop or nstop - 1...
   */

  double nu = 1e15; // TODO: set as 1e15 Hz for now.. update in future

  for (i = wind_start_index; i < wind_stop_index; i++)
  {
    cell_x = wmain[i].x;
    next_cell_x = wmain[i + 1].x;
    create_test_photon (&popac, nu, cell_x);

    ds = next_cell_x[0] - cell_x[0];
    opac_freq = radiation (&popac, ds); // / plasmamain[nplasma].density;
    Log ("Cell %i: ds %e opac_freq %e\n", i, ds, opac_freq);
  }

  free_opacity_arrays ();

  return EXIT_SUCCESS;
}

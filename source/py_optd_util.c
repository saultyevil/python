/** ************************************************************************* */
/**
 * @file     py_optical_depth_util.c
 * @author   Edward Parkinson
 * @date     May 2021
 *
 * @brief    Functions which aren't related to the transport of photons, or
 *           creation of the results.
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"
#include "py_optd.h"

/* ************************************************************************* */
/**
 * @brief  Initialize inclination angles for a spherical model
 *
 * @param[out]  SightLines_t *inclinations  The initialize inclinations structure
 *
 * @details
 *
 * Since this is for a 1D model, it does not matter which inclination angle
 * we use. We have opted to used an angle of 45 degrees, as this is also the
 * angle we use to define a 1D spherical grid. It should not matter either
 * way.
 *
 * ************************************************************************** */

static int
initialise_1d_angles (SightLines_t *inclinations)
{
  int len;
  const double phase = 1.0;
  const double default_angle = 45;
  CONFIG.n_inclinations = 1;

  inclinations[0].direction_vector[0] = sin (default_angle / RADIAN) * cos (-phase * 360.0 / RADIAN);
  inclinations[0].direction_vector[1] = sin (default_angle / RADIAN) * sin (-phase * 360.0 / RADIAN);
  inclinations[0].direction_vector[2] = cos (default_angle / RADIAN);

  len = snprintf (inclinations[0].name, NAMELEN, "A%02.0fP%04.2f", default_angle, phase);
  if (len < 0)
  {
    print_error ("there was an error writing the name to the sight lines array\n");
    exit (EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief  Initialize the inclination angles for a 2D outward run.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * @details
 *
 * If xxpsec has been read in, i.e. there were spectral cycles run the model,
 * then the same inclination angles in xxspec are used. Otherwise, a set of
 * default angles will be used instead, with a default phase of 0.5. The phase
 * should not matter, as it works for, i.e., 0.5 and 1 and gives the same
 * results as of writing this.
 *
 * ************************************************************************** */

static int
_mode_normal (struct SightLines *inclinations)
{
  int i;
  int len;
  const double phase = 1.0;
  const double default_angles[] = { 1.0, 10.0, 30.0, 45.0, 60.0, 75.0, 89.0 };
  const int n_default_angles = sizeof default_angles / sizeof default_angles[0];

  /*
   * Count the number of inclination angles provided by the -i command. Since
   * the array is initialized as -1, this assumes any values > -1 is a valid
   * inclination angle.
   */

  int n_user_input = 0;
  for (i = 0; i < MAX_CUSTOM_ANGLES; ++i)
  {
    if (CONFIG.inclinations[i] > -1)
    {
      ++n_user_input;
    }
  }

  // First of all use the user input
  if (n_user_input > 0)
  {
    CONFIG.n_inclinations = n_user_input;

    for (i = 0; i < n_user_input; i++)
    {
      len = snprintf (inclinations[i].name, NAMELEN, "A%02.0fP%04.2f", CONFIG.inclinations[i], phase);
      if (len < 0)
      {
        print_error ("there was an error writing the name to the sight lines array\n");
        exit (EXIT_FAILURE);
      }

      inclinations[i].direction_vector[0] = sin (CONFIG.inclinations[i] / RADIAN) * cos (-phase * 360.0 / RADIAN);
      inclinations[i].direction_vector[1] = sin (CONFIG.inclinations[i] / RADIAN) * sin (-phase * 360.0 / RADIAN);
      inclinations[i].direction_vector[2] = cos (CONFIG.inclinations[i] / RADIAN);
      inclinations[i].angle = CONFIG.inclinations[i];
    }
  }
  // Then use whatever is in the wind save and spec save
  else if (xxspec != NULL && geo.nangles > 0)
  {
    CONFIG.n_inclinations = geo.nangles;

    for (i = MSPEC; i < MSPEC + geo.nangles; i++)
    {
      strcpy (inclinations[i - MSPEC].name, xxspec[i].name);
      stuff_v (xxspec[i].lmn, inclinations[i - MSPEC].direction_vector);
      inclinations[i].angle = -1;       // todo: implement way to get angle xxspec
    }
  }
  // Otherwise, fall back to the default
  else
  {
    printf ("\nNo spec.save file has been found, using a default set of inclination angles\n\n");

    CONFIG.n_inclinations = n_default_angles;

    for (i = 0; i < n_default_angles; i++)
    {
      len = snprintf (inclinations[i].name, NAMELEN, "A%02.0fP%04.2f", default_angles[i], phase);
      if (len < 0)
      {
        print_error ("there was an error writing the name to the sight lines array\n");
        exit (EXIT_FAILURE);
      }

      inclinations[i].direction_vector[0] = sin (default_angles[i] / RADIAN) * cos (-phase * 360.0 / RADIAN);
      inclinations[i].direction_vector[1] = sin (default_angles[i] / RADIAN) * sin (-phase * 360.0 / RADIAN);
      inclinations[i].direction_vector[2] = cos (default_angles[i] / RADIAN);
      inclinations[i].angle = default_angles[i];
    }
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief  Initialize the inclination angles to find the photosphere for a
 *         2d model.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * @details
 *
 * The same function is called for both 1D and 2D models. This creates extra
 * work for 1D model, but as the algorithm takes very little time to run, it
 * does not matter.
 *
 * 500 inclination angles are defined, to very finely resolve the photosphere
 * surface. This is fixed for now. I think 500 is probably far too many, but
 * it takes absolutely no time to run. The results from 500 angles probably
 * need smoothing if the grid is coarse.
 *
 * ************************************************************************** */

static int
_mode_surface (struct SightLines *inclinations)
{
  const double phase = 1.0;
  const double d_theta = 90.0 / (float) MAX_ANGLES;

  CONFIG.n_inclinations = MAX_ANGLES;

  for (int i = 0; i < MAX_ANGLES; i++)
  {
    double angle = i * d_theta;
    inclinations[i].direction_vector[0] = sin (angle / RADIAN) * cos (-phase * 360.0 / RADIAN);
    inclinations[i].direction_vector[1] = sin (angle / RADIAN) * sin (-phase * 360.0 / RADIAN);
    inclinations[i].direction_vector[2] = cos (angle / RADIAN);
    inclinations[i].angle = angle;
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

static void
initialise_2d_angles (struct SightLines *inclinations)
{
  switch (CONFIG.run_mode)
  {
  case MODE_SURFACE:
    _mode_surface (inclinations);
    break;
  default:
    _mode_normal (inclinations);
    break;
  }
}

/* ************************************************************************* */
/**
 * @brief  Wrapper function for initializing the inclination angles depending
 *         on the run mode of the program.
 *
 * @param[out]  int n_angles  The number of angles initialized
 *
 * @return  SightLines_t *inclinations  The initialize inclinations structure
 *
 * ************************************************************************** */

int
initialize_inclination_angles (struct SightLines *inclinations)
{
  // inclinations need to be malloc'd already, I think.

  if (zdom[CONFIG.domain].coord_type == SPHERICAL)
  {
    initialise_1d_angles (inclinations);
  }
  else
  {
    initialise_2d_angles (inclinations);
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief  Generate a photon packet with a given frequency nu for use with the
 *         optical depth diagnostic routines.
 *
 * @param[in,out]  p_out  The photon packet to initialise
 * @param[in]  freq  The frequency of the photon packet
 *
 * @return  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * This routine assumes that the user wishes for the photon to be generated from
 * the surface of the central source, taking into account the current direction
 * of the sight-line being extracted. It's currently not possible, without a
 * tedious re-write, to place the photon at the very origin of the grid when
 * there is a central source because of how we check boundaries.
 *
 * Note that photons are initialised with a weight of f_tot as photons are
 * required to have weight, but since functions do not care about the weight of
 * the photon, it is set to something large to make sure it does not get
 * destroyed by accident somehow.
 *
 * ************************************************************************** */

int
initialise_photon_packet (PhotPtr photon, double frequency, double *direction)
{
  static int np = 0;

  if (frequency < 0)
  {
    print_error ("photon can't be created with negative frequency\n");
    exit (EXIT_FAILURE);
  }

  photon->np = ++np;
  photon->freq = photon->freq_orig = frequency;
  photon->origin = photon->origin_orig = PTYPE_STAR;
  photon->istat = P_INWIND;
  photon->w = photon->w_orig = geo.f_tot;
  photon->tau = 0.0;
  photon->frame = F_OBSERVER;
  photon->x[0] = photon->x[1] = photon->x[2] = 0.0;
  stuff_v (direction, photon->lmn);

  switch (CONFIG.run_mode)
  {
  case MODE_SURFACE:           // Move to edge of wind and point it inward
    move_phot (photon, zdom[CONFIG.domain].rmax - DFUDGE);
    for (int i = 0; i < 3; ++i)
    {
      photon->lmn[i] *= -1.0;
    }
    break;
  case MODE_CELL_SPECTRUM:     // Move to the bottom left of the cell
    photon->x[0] = 56000000000000;
    photon->x[2] = 11300000000000;
    break;
  default:                     // Move to the inner edge of the outflow
    move_phot (photon, geo.rstar + DFUDGE);
    break;
  }

  photon->grid = where_in_grid (CONFIG.domain, photon->x);

  return EXIT_SUCCESS;
}

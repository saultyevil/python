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
#include <unistd.h>

#include "atomic.h"
#include "python.h"
#include "py_optd.h"

/* ************************************************************************** */
/**
 * @brief Initialize a Sirocco/Python wind structure and other associated data.
 *
 * @details
 *
 * Data is read in from the wind_save and spec_save files. This function uses
 * some of Sirocco's global variables and functions to read in the data. This is
 * why there are no arguments to this function, as Sirocco handles all of this
 * with global state.
 *
 * ************************************************************************** */

int
initialize_wind_structures (void)
{
  char windsave_filename[LINELENGTH + LINELENGTH];
  char specsave_filename[LINELENGTH + LINELENGTH];
  snprintf (windsave_filename, LINELENGTH + LINELENGTH, "%s.wind_save", files.root);
  snprintf (specsave_filename, LINELENGTH + LINELENGTH, "%s.spec_save", files.root);

  /*
   * Read in the wind_save file and initialize the wind cones and DFUDGE which
   * are important for photon transport. The atomic data is also read in at
   * this point (which is also very important)
   */

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    print_error ("Failed to allocate memory for wind domain(s)\n");
    return EXIT_FAILURE;
  }

  if (wind_read (windsave_filename) < 0)
  {
    print_error ("Failed to read in wind save file from %s\n", windsave_filename);
    return EXIT_FAILURE;
  }

  DFUDGE = setup_dfudge ();
  setup_windcone ();

  /*
   * If a spec_save exists, and there are spectral cycles (possibly a redundant
   * check), then read in the spec_save file.
   */

  if (access (specsave_filename, F_OK) == 0)
  {
    if (geo.pcycle > 0)
    {
      if (spec_read (specsave_filename) < 0)
      {
        print_error ("Failed to open %s, when it should exist for this simulation\n", specsave_filename);
        return EXIT_FAILURE;
      }
    }
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************** */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

void
set_frequency_range (void)
{
  double freq_min = 0;
  double freq_max = 0;
  double cli_freq_min = CONFIG.arg_freq_min;
  double cli_freq_max = CONFIG.arg_freq_max;

  if (cli_freq_min > 0 || cli_freq_max > 0)
  {
    if (cli_freq_min > 0)
    {
      freq_min = cli_freq_min;
    }
    else
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
    }

    if (cli_freq_max > 0)
    {
      freq_max = cli_freq_max;
    }
    else
    {
      freq_max = VLIGHT / (100 * ANGSTROM);
    }

    if (freq_max < freq_min)
    {
      print_error ("frequency range given has set freq_max (%e) < freq_min (%e) \n", freq_max, freq_min);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    if ((geo.nangles == 0 && xxspec == NULL) || (geo.swavemax == 0 && geo.swavemin == 0))
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
      freq_max = VLIGHT / (100 * ANGSTROM);
    }
    else
    {
      freq_min = VLIGHT / (geo.swavemax * ANGSTROM);
      freq_max = VLIGHT / (geo.swavemin * ANGSTROM);
      if (sane_check (freq_min))
      {
        freq_min = VLIGHT / (10000 * ANGSTROM);
        print_error ("freq_min has an invalid value setting to %e\n", freq_min);
      }
      if (sane_check (freq_max))
      {
        freq_max = VLIGHT / (100 * ANGSTROM);
        print_error ("freq_min has an invalid value setting to %e\n", freq_max);
      }
    }
  }

  CONFIG.arg_freq_min = freq_min;
  CONFIG.arg_freq_max = freq_max;
}

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
  CONFIG.arg_num_inc = 1;

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
initialise_2d_normal (struct SightLines *inclinations)
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
    if (CONFIG.arg_inclinations[i] > -1)
    {
      ++n_user_input;
    }
  }

  // First of all use the user input
  if (n_user_input > 0)
  {
    CONFIG.arg_num_inc = n_user_input;

    for (i = 0; i < n_user_input; i++)
    {
      len = snprintf (inclinations[i].name, NAMELEN, "A%02.0fP%04.2f", CONFIG.arg_inclinations[i], phase);
      if (len < 0)
      {
        print_error ("there was an error writing the name to the sight lines array\n");
        exit (EXIT_FAILURE);
      }

      inclinations[i].direction_vector[0] = sin (CONFIG.arg_inclinations[i] / RADIAN) * cos (-phase * 360.0 / RADIAN);
      inclinations[i].direction_vector[1] = sin (CONFIG.arg_inclinations[i] / RADIAN) * sin (-phase * 360.0 / RADIAN);
      inclinations[i].direction_vector[2] = cos (CONFIG.arg_inclinations[i] / RADIAN);
      inclinations[i].angle = CONFIG.arg_inclinations[i];
    }
  }
  // Then use whatever is in the wind save and spec save
  else if (xxspec != NULL && geo.nangles > 0)
  {
    CONFIG.arg_num_inc = geo.nangles;

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

    CONFIG.arg_num_inc = n_default_angles;

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
initialise_2d_surface (struct SightLines *inclinations)
{
  const double phase = 1.0;
  const double d_theta = 90.0 / (float) MAX_ANGLES;

  CONFIG.arg_num_inc = MAX_ANGLES;

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
  switch (CONFIG.mode)
  {
  case MODE_FIND_SURFACE:
    initialise_2d_surface (inclinations);
    break;
  default:
    initialise_2d_normal (inclinations);
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

  switch (CONFIG.mode)
  {
  case MODE_FIND_SURFACE:      // Move to edge of wind and point it inward
    move_phot (photon, zdom[CONFIG.domain].rmax - DFUDGE);
    for (int i = 0; i < 3; ++i)
    {
      photon->lmn[i] *= -1.0;
    }
    break;
  case MODE_CELL_SPECTRUM:     // Move to the bottom left of the cell -- this is not great, passing data with global state!
    photon->x[0] = wmain[CONFIG.arg_wind_elem].x[0];
    photon->x[1] = wmain[CONFIG.arg_wind_elem].x[1];
    photon->x[2] = wmain[CONFIG.arg_wind_elem].x[2];
    break;
  default:                     // Move to the inner edge of the outflow
    move_phot (photon, geo.rstar + DFUDGE);
    break;
  }

  int wind_status = where_in_wind (photon->x, &CONFIG.domain);
  if (wind_status < 0)
  {
    return EXIT_FAILURE;
  }
  photon->grid = where_in_grid (CONFIG.domain, photon->x);

  return EXIT_SUCCESS;
}

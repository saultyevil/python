/** ************************************************************************* */
/**
 * @file     py_optical_depth_transport.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    Functions related to the transport of photons and the main
 *           algorithm functions.
 *
 * ************************************************************************** */

#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"
#include "py_optd.h"

#define MAX_DIFF (VCHECK / VLIGHT)      // For linear velocity requirement for photon transport
#define MAX_TRANSLATE_IN_SPACE 10

/* ************************************************************************* */
/**
 * @brief  Move a photon across a cell, whilst calculating the optical depth
 *         and column density
 *
 * @param[in] photon  The photon packet to move
 * @param[in,out] *column_density_out  The column density moved through, updated
 *                                     in place
 * @param[in,out] *optical_depth_out   The optical depth moved through, updated
 *                                     in place
 *
 * @return p_istat  The current photon status or EXIT_FAILURE on failure.
 *
 * @details
 *
 * A photon is moved some distance along a cell in its current direction. This
 * may be the entire cell, but more likely it is some smaller path which ensures
 * the frequency shift is not too large that a linear approximation of the
 * frequency change is invalid.
 *
 * ************************************************************************** */

int
move_photon_across_cell (PhotPtr photon, double *column_density_out, double *optical_depth_out)
{
  int photon_status;
  double kappa_total, density;
  double smax, diff;
  double freq_inner, freq_outer, mean_freq;
  struct photon photon_start, photon_stop, photon_now;

  int n_dom;
  int n_plasma;
  WindPtr wind_cell;
  PlasmaPtr plasma_cell;

  wind_cell = &wmain[photon->grid];
  n_dom = wind_cell->ndom;
  n_plasma = wmain[photon->grid].nplasma;
  plasma_cell = &plasmamain[n_plasma];

  if (CONFIG.column_density == COLUMN_MODE_RHO)
  {
    density = plasma_cell->rho;
  }
  else
  {
    density = plasma_cell->density[CONFIG.column_density_ion_number];
  }

  // We use SMAX_FRAC to prevent photons moving too far which can, in some
  // cases, create weird looking edges or miss where a optical depth surface is
  // reached.
  smax = smax_in_cell (photon) * SMAX_FRAC;
  if (smax < -DBL_EPSILON)
  {
    print_error ("smax %e < 0 in cell %d\n", smax, photon->grid);
    return EXIT_FAILURE;
  }

  // Transform the photon into the CMF and create a photon in the CMF at the
  // end of the path of length smax
  observer_to_local_frame (photon, &photon_start);
  stuff_phot (photon, &photon_stop);
  move_phot (&photon_stop, smax);
  observer_to_local_frame (&photon_stop, &photon_stop);

  // At this point p_start and p_stop are in the local frame
  // at the and p_stop is at the maximum distance it can
  // travel. We want to check that the frequency shift is
  // not too great along the path that a linear approximation
  // to the change in frequency is not reasonable

  while (smax > DFUDGE)
  {
    stuff_phot (photon, &photon_now);
    move_phot (&photon_now, smax * 0.5);
    observer_to_local_frame (&photon_now, &photon_now);
    diff = fabs (photon_now.freq - 0.5 * (photon_start.freq + photon_stop.freq)) / photon_start.freq;
    if (diff < MAX_DIFF)
      break;
    stuff_phot (&photon_now, &photon_stop);
    smax *= 0.5;
  }

  freq_inner = photon_start.freq;
  freq_outer = photon_stop.freq;
  mean_freq = 0.5 * (freq_inner + freq_outer);

  // Now we can finally calculate the opacity due to all the continuum
  // processes. In macro-atom mode, we need to calculate the continuum opacity
  // using kappa_bf and kappa_ff using the macro treatment. For simple mode, we
  // can simply use radiation which **SHOULD** return the continuum opacity
  // as well, plus something from induced Compton heating. In either cases,
  // we still then need to add the optical depth from electron scattering at
  // the end.

  kappa_total = 0;

  if (CONFIG.mode != MODE_SURFACE)
  {
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      kappa_total += radiation (photon, smax);
    }
    else                        // macro atom case
    {
      if (wind_cell->vol > 0)
      {
        kappa_total += kappa_bf (plasma_cell, freq_inner, 0);
        kappa_total += kappa_ff (plasma_cell, freq_inner);
      }
    }
  }
  if (CONFIG.ignore_electron_scattering == false)
  {
    kappa_total += klein_nishina (mean_freq) * plasma_cell->ne * zdom[n_dom].fill;
  }

  *column_density_out += smax * density;
  *optical_depth_out += smax * kappa_total;
  move_phot (photon, smax);
  photon_status = photon->istat;

  return photon_status;
}

/* ************************************************************************* */
/**
 * @brief Calculate the optical depth across a cell for the given photon
 *
 * @param[in]  photon_in
 * @param[out] *optical_depth_out
 *
 *
 * @return  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * ************************************************************************** */

int
integrate_tau_across_cell (PhotPtr photon, double *column_density_out, double *optical_depth_out)
{
  int ndom;

  where_in_wind (photon->x, &ndom);
  int grid_start = where_in_grid (ndom, photon->x);

  double column_density = 0;
  double optical_depth = 0;

  while (photon->grid == grid_start)
  {
    if (photon->grid < 0 || photon->grid > zdom[ndom].ndim2 - 1)
    {
      print_error ("Photon is not in the grid, aborting this photon\n");
      return EXIT_FAILURE;
    }
    // move across cell and calculate quantities at the same time
    int error = move_photon_across_cell (photon, &column_density, &optical_depth);
    if (error)
    {
      return EXIT_FAILURE;
    }
    // check to see if it's still in the same cell
    photon->grid = where_in_grid (ndom, photon->x);
  }

  *column_density_out = column_density;
  *optical_depth_out = optical_depth;

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief           Extract the optical depth the photon packet porig must
 *                  travel through to reach the observer.
 *
 * @param[in]  photon  The photon packet to extract
 * @param[out]  *c_column_density  The column depth of the extracted photon angle
 * @param[out]  *c_optical_depth  The optical depth from photon origin to
 *                                the observer
 *
 * @return  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * The direction of the observer is set in porig.lmn and should be set prior
 * to passing a photon to this function. If any values of lmn are negative, the
 * photon will translate in the negative direction. However, have no fear as this
 * is normal and is fine due to the assumed symmetry of models in Python.
 *
 * ************************************************************************** */

int
integrate_tau_across_wind (PhotPtr photon_in, double *column_density, double *optical_depth)
{
  int ndom;
  double _norm[3];
  struct photon photon;

  enum istat_enum photon_status = P_INWIND;
  stuff_phot (photon_in, &photon);      // we need a copy of the original for walls
  int num_times_in_space = 0;

  while (photon_status == P_INWIND)
  {
    int wind_status = where_in_wind (photon.x, &ndom);

    if (wind_status < 0)        // this shouldn't happen, so we'll limit the number of times we translate in space
    {
      num_times_in_space += 1;
      translate_in_space (&photon);
      if (num_times_in_space > MAX_TRANSLATE_IN_SPACE)
      {
        print_error ("Photon has been outside the grid too many times, something is wrong");
        return EXIT_FAILURE;
      }
    }
    else
    {
      int error = integrate_tau_across_cell (&photon, column_density, optical_depth);
      if (error)
      {
        return EXIT_FAILURE;
      }
    }

    if (CONFIG.mode == MODE_SURFACE)
    {
      if (photon.tau > CONFIG.arg_tau_surface)
      {
        photon_status = P_ABSORB;
        break;
      }
    }

    photon_status = walls (&photon, photon_in, _norm);
  }

  // If a photon hits the surface, something will have gone wrong. Throw that
  // photon away and return an error
  if (CONFIG.mode == MODE_SPECTRUM)
  {
    if (photon_status == P_HIT_STAR || photon_status == P_HIT_DISK)
    {
      print_error ("photon hit central source or disk surface, when it shouldn't be able to\n");
      return EXIT_FAILURE;
    }
  }
  if (CONFIG.mode == MODE_SURFACE)
  {
    if (photon_status == P_HIT_DISK)
    {
      print_error ("photon hit disk surface, when it should hit the central source\n");
      return EXIT_FAILURE;
    }
  }

  // Some modes will re-use the photon to do some additional stuff
  stuff_phot (&photon, photon_in);

  return EXIT_SUCCESS;
}

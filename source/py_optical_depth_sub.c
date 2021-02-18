/** ************************************************************************* */
/**
 * @file     py_optical_depth_sub.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    The main resting place for most functions related to providing
 *           optical depth diagnostics.
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"
#include "py_optical_depth.h"

const PIEdges_t PI_EDGES_TO_MEASURE[] = {
  {"HLymanEdge", 3.387485e+15},
  {"HBalmerEdge", 8.293014e+14},
  {"HeI24eVEdge", 5.9483e+15},
  {"HeII54eVEdge", 1.394384e+16}
};

/* ************************************************************************** */
/**
 * @brief  Print the various optical depths calculated using this routine
 *
 * @param[in]  optical_depths  The 2d array containing the optical depth for
 *                             each observer and tau
 * @param[in]  column_densities  The 2d array containing the column densities
 *                               for each observer
 *
 * @details
 *
 * Prints the different optical depths for each angle and optical depth.
 * Historically, this used to also print to file. But I found that it was mostly
 * useless to do this.
 *
 * ************************************************************************** */

void
print_optical_depths (const double *optical_depth_values, const double *column_density_values)
{
  int i, j;
  int c_linelen;
  char str[LINELENGTH];
  const int MAX_COL = 120;

  printf ("\nOptical depths along the defined line of sights for domain %i:\n\n", DOMAIN_TO_CONSIDER);

  for (i = 0; i < N_INCLINATION_ANGLES; i++)
  {
    if (COLUMN_MODE == COLUMN_MODE_RHO)
    {
      printf ("%-8s: Mass column density     : %3.2e cm^-2\n", INCLINATION_ANGLES[i].name, column_density_values[i]);
      printf ("%-8s: Hydrogen column density : %3.2e cm^-2\n", "", column_density_values[i] * rho2nh);
    }
    else
    {
      printf ("%-8s: %s %i column density    : %3.2e cm^-2\n", INCLINATION_ANGLES[i].name, ele[ion[COLUMN_MODE_ION_NUMBER].nelem].name,
              ion[COLUMN_MODE_ION_NUMBER].istate, column_density_values[i]);
    }

    c_linelen = 0;
    for (j = 0; j < N_PI_EDGES; j++)
    {
      c_linelen += sprintf (str, "tau_%-9s: %3.2e  ", PI_EDGES_TO_MEASURE[j].name, optical_depth_values[i * N_PI_EDGES + j]);
      if (c_linelen > MAX_COL)
      {
        c_linelen = 0;
        printf ("\n");
      }
      printf ("%s", str);
    }
    printf ("\n\n");
  }
}

/* ************************************************************************* */
/**
 * @brief           Write the various optical depth spectra to file
 *
 * @param[in]  tau_spectrum  The various optical depth spectra
 * @param[in]  freq_min  The smallest frequency in the spectra
 * @param[in]  d_freq  The size of the frequncy bins
 *
 * @details
 *
 * Simply write the optical depth spectra to the file named root.tau_spec.diag.
 * This file will be located in the diag folder.
 *
 * ************************************************************************** */

void
write_optical_depth_spectrum (const double *tau_spectrum, const double freq_min, const double d_freq)
{
  int i, j;
  double c_wavelength, c_frequency;
  char filename[LINELENGTH + 12];
  FILE *fp;

  sprintf (filename, "%s.spec_tau", files.root);

  if ((fp = fopen (filename, "w")) == NULL)
  {
    printf ("write_optical_depth_spectrum: uh oh, could not open optical depth spectrum output file\n");
    exit (EXIT_FAILURE);
  }

  /*
   * Write out the tau spectrum header
   */

  fprintf (fp, "%-12s %-12s ", "Freq.", "Lambda");
  for (j = 0; j < N_INCLINATION_ANGLES; j++)
    fprintf (fp, "%-12s ", INCLINATION_ANGLES[j].name);
  fprintf (fp, "\n");

  /*
   * Write out the tau spectrum for each inclination angle
   */

  c_frequency = freq_min;
  for (i = 0; i < N_FREQ_BINS; i++)
  {
    c_wavelength = VLIGHT / c_frequency / ANGSTROM;
    fprintf (fp, "%-12e %-12e ", c_frequency, c_wavelength);
    for (j = 0; j < N_INCLINATION_ANGLES; j++)
      fprintf (fp, "%-12e ", tau_spectrum[j * N_FREQ_BINS + i]);
    fprintf (fp, "\n");
    c_frequency += d_freq;
  }

  fflush (fp);                  // probably not required, but was needed when this was part of Python rather than standalone

  if (fclose (fp))
    printf ("write_optical_depth_spectrum: uh oh, could not close optical depth spectrum output file\n");
}

/* ************************************************************************* */
/**
 * @brief  Calculate the total optical depth a photon experiences across the
 *         cell of distance SMAX_FRAC * smax.
 * @param[in]  photon  The photon packet
 * @param[in,out]  *c_column_density  The column density the photon has moved
 *                                    through
 * @param[in,out]  *c_optical_depth  The optical depth experienced by the photon
 *
 * @return p_istat  The current photon status or EXIT_FAILURE on failure.
 *
 * @details
 *
 * This function is concerned with finding the opacity of the photon's current
 * cell, the distance the photon can move in the cell and hence it increments
 * the optical depth tau a photon has experienced as it moves through the wind.
 *
 * ************************************************************************** */

const double MAXDIFF = VCHECK / VLIGHT; // For linear velocity requirement for photon transport

int
calculate_tau_across_cell (PhotPtr photon, double *c_column_density, double *c_optical_depth)
{
  int p_istat;
  int n_dom, n_plasma;
  double kappa_total, density;
  double smax;
  double diff;
  double freq_inner, freq_outer, mean_freq;
  WindPtr c_wind_cell;
  PlasmaPtr c_plasma_cell;
  struct photon p_start, p_stop, p_now;

  if ((photon->grid = where_in_grid (wmain[photon->grid].ndom, photon->x)) < 0)
  {
    printf ("calculate_tau_across_cell: photon is not in grid in cell %i\n", photon->grid);
    return EXIT_FAILURE;
  }

  /*
   * Determine where in the plasma grid the cell is. This is required so a
   * column density can be calculated. Also calculate the maximum distance the
   * photon can traverse across the cell
   */

  c_wind_cell = &wmain[photon->grid];
  n_dom = c_wind_cell->ndom;
  n_plasma = wmain[photon->grid].nplasma;
  c_plasma_cell = &plasmamain[n_plasma];

  if (COLUMN_MODE == COLUMN_MODE_RHO)
  {
    density = c_plasma_cell->rho;
  }
  else
  {
    density = c_plasma_cell->density[COLUMN_MODE_ION_NUMBER];
  }

  smax = smax_in_cell (photon);
  if (smax < 0)
  {
    printf ("calculate_tau_across_cell: abnormal value of smax for photon\n");
    return EXIT_FAILURE;
  }

  /*
   * Transform the photon into the CMF and create a photon in the CMF at the
   * end of the path of length smax
   */

  observer_to_local_frame (photon, &p_start);
  stuff_phot (photon, &p_stop);
  move_phot (&p_stop, smax);
  observer_to_local_frame (&p_stop, &p_stop);

  /* At this point p_start and p_stop are in the local frame
   * at the and p_stop is at the maximum distance it can
   * travel. We want to check that the frequency shift is
   * not too great along the path that a linear approximation
   * to the change in frequency is not reasonable
   */

  while (smax > DFUDGE)
  {
    stuff_phot (photon, &p_now);
    move_phot (&p_now, smax * 0.5);
    observer_to_local_frame (&p_now, &p_now);
    diff = fabs (p_now.freq - 0.5 * (p_start.freq + p_stop.freq)) / p_start.freq;
    if (diff < MAXDIFF)
      break;
    stuff_phot (&p_now, &p_stop);
    smax *= 0.5;
  }

  freq_inner = p_start.freq;
  freq_outer = p_stop.freq;
  mean_freq = 0.5 * (freq_inner + freq_outer);

  /*
   * Now we can finally calculate the opacity due to all the continuum
   * processes. In macro-atom mode, we need to calculate the continuum opacity
   * using kappa_bf and kappa_ff using the macro treatment. For simple mode, we
   * can simply use radiation which **SHOULD** return the continuum opacity
   * as well, plus something from induced Compton heating. In either cases,
   * we still then need to add the optical depth from electron scattering at
   * the end.
   */

  kappa_total = 0;

  if (geo.rt_mode == RT_MODE_2LEVEL)
  {
    radiation (photon, smax, &kappa_total);
  }
  else                          // macro atom case
  {
    if (c_wind_cell->vol > 0)
    {
      kappa_total += kappa_bf (c_plasma_cell, freq_inner, 0);
      kappa_total += kappa_ff (c_plasma_cell, freq_inner);
    }
  }

  kappa_total += klein_nishina (mean_freq) * c_plasma_cell->ne * zdom[n_dom].fill;

  /*
   * Increment the optical depth and column density variables and move the
   * photon to the edge of the cell
   */

  *c_column_density += smax * density;
  *c_optical_depth += smax * kappa_total;
  move_phot (photon, smax);
  p_istat = photon->istat;

  return p_istat;
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
extract_tau (PhotPtr photon, double *c_column_density, double *c_optical_depth)
{
  int n_dom;
  enum istat_enum p_istat;
  const int max_translate_in_space = 10;
  int n_in_space;
  double norm[3];
  struct photon p_extract;

  p_istat = P_INWIND;           // assume photon is in wind for initialisation reasons
  stuff_phot (photon, &p_extract);

  n_in_space = 0;
  while (p_istat == P_INWIND)
  {
    if (where_in_wind (p_extract.x, &n_dom) < 0)
    {
      translate_in_space (&p_extract);
      if (++n_in_space > max_translate_in_space)
      {
        printf ("extract_tau: tau_extract photon transport ended due to too many translate_in_space\n");
        return EXIT_FAILURE;
      }
    }
    else if ((p_extract.grid = where_in_grid (n_dom, p_extract.x)) >= 0)
    {
      if (calculate_tau_across_cell (&p_extract, c_column_density, c_optical_depth))
        return EXIT_FAILURE;
    }
    else
    {
      printf ("extract_tau: photon in unknown location grid stat %i\n", p_extract.grid);
      return EXIT_FAILURE;
    }

    p_istat = walls (&p_extract, photon, norm);
  }

  if (p_istat == P_HIT_STAR || p_istat == P_HIT_DISK)
  {
    printf ("extract_tau: photon hit central source or disk incorrectly istat = %i\n", p_istat);
    return EXIT_FAILURE;
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
create_photon (PhotPtr p_out, double freq, double *lmn)
{
  if (freq < 0)
  {
    printf ("create_photon: photon can't be created with negative frequency\n");
    return EXIT_FAILURE;
  }

  stuff_v (lmn, p_out->lmn);
  p_out->freq = p_out->freq_orig = freq;
  p_out->origin = p_out->origin_orig = PTYPE_DISK;
  p_out->istat = P_INWIND;
  p_out->w = p_out->w_orig = geo.f_tot;
  p_out->tau = 0.0;
  p_out->x[0] = p_out->x[1] = p_out->x[2] = 0;
  p_out->frame = F_OBSERVER;
  move_phot (p_out, geo.rstar + DFUDGE);

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief  Initialise the viewing angles for the optical depth diagnostic
 *         routines.
 *
 * @details
 *
 * This purpose of this function is to initialise the angles of which to extract
 * optical depth information from. If xxspec is initialised, then the optical
 * depth information will be extracted from the observer angles provided in the
 * parameter file for spectrum generation. Otherwise, the routines will instead
 * use a set of pre-defined angles.
 *
 * ************************************************************************** */

void
initialize_inclination_angles (void)
{
  int i;
  long mem_req;
  const double default_phase = 0.5;
  const double default_angles[] = { 0.0, 10.0, 30.0, 45.0, 60.0, 75.0, 85.0, 90.0 };
  const int n_default_angles = sizeof default_angles / sizeof default_angles[0];

  /*
   * Use the angles specified for by the user for spectrum generation, this
   * requires for xxspec to be initialised.
   */

  if (geo.nangles > 0 && xxspec != NULL)
  {
    N_INCLINATION_ANGLES = geo.nangles;
    INCLINATION_ANGLES = calloc (geo.nangles, sizeof *INCLINATION_ANGLES);
    if (INCLINATION_ANGLES == NULL)
    {
      mem_req = geo.nangles * (int) sizeof *INCLINATION_ANGLES;
      printf ("initialize_inclination_angles: cannot allocate %ld bytes for observers array\n", mem_req);
      exit (EXIT_FAILURE);
    }
    else                        // Use else to avoid compiler warning
    {
      for (i = MSPEC; i < MSPEC + geo.nangles; i++)
      {
        strcpy (INCLINATION_ANGLES[i - MSPEC].name, xxspec[i].name);
        stuff_v (xxspec[i].lmn, INCLINATION_ANGLES[i - MSPEC].lmn);
      }
    }
  }

  /*
   * If no spectrum cycles are being run, or if the model is being restarted
   * without initialising xxspec, then a default set of angles is used instead.
   */

  else
  {
    printf ("tau_spectrum: as there are no spectrum cycles or observers defined, a set of default angles will be used instead\n");
    N_INCLINATION_ANGLES = n_default_angles;
    INCLINATION_ANGLES = calloc (n_default_angles, sizeof *INCLINATION_ANGLES);
    if (INCLINATION_ANGLES == NULL)
    {
      mem_req = n_default_angles * (int) sizeof *INCLINATION_ANGLES;
      printf ("initialize_inclination_angles: cannot allocate %ld bytes for observers array\n", mem_req);
      exit (EXIT_FAILURE);
    }
    else
    {
      for (i = 0; i < n_default_angles; i++)
      {
        sprintf (INCLINATION_ANGLES[i].name, "A%02.0fP%04.2f", default_angles[i], default_phase);
        INCLINATION_ANGLES[i].lmn[0] = sin (default_angles[i] / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
        INCLINATION_ANGLES[i].lmn[1] = sin (default_angles[i] / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
        INCLINATION_ANGLES[i].lmn[2] = cos (default_angles[i] / RADIAN);
      }
    }
  }
}

/* ************************************************************************* */
/**
 * @brief  Create spectra of tau vs lambda for each observer angle
 *
 * @details
 *
 * This is the main function which will generate the optical depth spectra for
 * each observer angle in xxspec. The algorithm is similar to extract and the
 * tau_diag algorithm which this function is called in.
 *
 * A photon is generated at the central source of the model and is extracted
 * from this location towards the observer where it escapes, where extract_tau
 * returns the integrated optical depth along its path to escape. This is done
 * for a range of photon frequencies to find how optical depth changes with
 * frequency.
 *
 * This processes can take some time compared to tau_evalulate_photo_edges. But,
 * since N_FREQ_BINS photons are being generated for each spectrum and the fact
 * that these photons do not interact, the spectra does not actually take that
 * long to complete.
 *
 * ************************************************************************** */

void
create_optical_depth_spectrum (void)
{
  int i, j;
  double *tau_spectrum;
  double c_optical_depth, c_column_density;
  double c_frequency, freq_min, freq_max, d_freq;
  struct photon photon;

  printf ("Creating optical depth spectra:\n");

  if ((tau_spectrum = calloc (N_INCLINATION_ANGLES * N_FREQ_BINS, sizeof *tau_spectrum)) == NULL)
  {
    printf ("create_optical_depth_spectrum: cannot allocate %lu bytes for tau_spectrum\n",
            N_INCLINATION_ANGLES * N_FREQ_BINS * sizeof *tau_spectrum);
    return;
  }

  /*
   * Define the limits of the spectra in frequency space. If xxpsec is NULL,
   * then the frequency range will be over a default 100 - 10,000 Angstrom
   * band.
   */

  if ((geo.nangles == 0 && xxspec == NULL) || (geo.swavemax == 0 && geo.swavemin == 0))
  {
    printf ("create_optical_depth_spectrum: xxspec is uninitialized, defaulting spectral wavelength range to 100 - 10,000 Angstrom\n");
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
      printf ("create_optical_depth_spectrum: freq_min has an invalid value setting to %e\n", freq_min);
    }
    if (sane_check (freq_max))
    {
      freq_max = VLIGHT / (100 * ANGSTROM);
      printf ("create_optical_depth_spectrum: freq_min has an invalid value setting to %e\n", freq_max);
    }
  }

  d_freq = (freq_max - freq_min) / N_FREQ_BINS;
  kbf_need (freq_min, freq_max);

  /*
   * Now create the optical depth spectra for each inclination
   */

  for (i = 0; i < N_INCLINATION_ANGLES; i++)
  {
    printf ("  - Creating spectrum: %s\n", INCLINATION_ANGLES[i].name);
    c_frequency = freq_min;
    for (j = 0; j < N_FREQ_BINS; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;
      if ((create_photon (&photon, c_frequency, INCLINATION_ANGLES[i].lmn)) == EXIT_FAILURE)
      {
        printf ("create_optical_depth_spectrum: skipping photon of frequency %e when creating spectrum\n", c_frequency);
        continue;
      }
      if ((extract_tau (&photon, &c_column_density, &c_optical_depth)) == EXIT_FAILURE)
        continue;
      tau_spectrum[i * N_FREQ_BINS + j] = c_optical_depth;
      c_frequency += d_freq;
    }
  }

  write_optical_depth_spectrum (tau_spectrum, freq_min, d_freq);
  if (tau_spectrum)             // I didn't think I needed to do this, but if I don't, there are runtime errors
    free (tau_spectrum);
}

/* ************************************************************************* */
/**
 * @brief Calculate the optical depth for various optical depth edges and
 *        extract the column density.
 *
 * @details
 *
 * This is the main function which will control the procedure for calculating
 * various diagnostic numbers for the optical depth's experienced in the current
 * model. Namely, this function aims to show the total integrated optical depth
 * to each observer angle using (originally) the following optical depths:
 *
 *  - Lymann edge
 *  - Balmer edge
 *  - Helium II edge
 *
 * Once these integrated optical depths have been calculated for each angle, a
 * spectrum of optical depth vs wavelength is created for each angle.
 *
 * The aim of these diagnostic numbers it to provide some sort of quick metric
 * on the optical thickness of the current model.
 *
 * ************************************************************************** */

void
calculate_PI_edge_optical_depth (void)
{
  int i, j;
  double c_frequency, c_optical_depth, c_column_density;
  double *optical_depth_values, *column_density_values;
  struct photon photon;

  if ((optical_depth_values = calloc (N_INCLINATION_ANGLES * N_PI_EDGES, sizeof *optical_depth_values)) == NULL)
  {
    printf ("calculate_PI_edge_optical_depth: cannot allocate %lu bytes for optical_depths\n",
            N_INCLINATION_ANGLES * N_PI_EDGES * sizeof *optical_depth_values);
    return;
  }
  if ((column_density_values = calloc (N_INCLINATION_ANGLES, sizeof *optical_depth_values)) == NULL)
  {
    printf ("calculate_PI_edge_optical_depth: cannot allocate %lu bytes for column_densities\n",
            N_INCLINATION_ANGLES * sizeof *optical_depth_values);
    free (optical_depth_values);
    return;
  }

  /*
   * Now extract the optical depths and mass column densities. We loop over
   * each PI edge for each inclination angle.
   */

  for (i = 0; i < N_INCLINATION_ANGLES; i++)
  {
    for (j = 0; j < N_PI_EDGES; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;
      c_frequency = PI_EDGES_TO_MEASURE[j].freq;

      if (create_photon (&photon, c_frequency, INCLINATION_ANGLES[i].lmn) == EXIT_FAILURE)
      {
        printf ("calculate_PI_edge_optical_depth: skipping photon of frequency %e\n", c_frequency);
        continue;
      }

      if (extract_tau (&photon, &c_column_density, &c_optical_depth) == EXIT_FAILURE)
        continue;               // do not throw extra warning when one is already thrown in extract_tau

      optical_depth_values[i * N_PI_EDGES + j] = c_optical_depth;
      column_density_values[i] = c_column_density;
    }
  }

  print_optical_depths (optical_depth_values, column_density_values);

  /*
   * The if's here are to avoid warnings with potentially freeing already
   * deallocated memory
   */

  if (optical_depth_values)
    free (optical_depth_values);
  if (column_density_values)
    free (column_density_values);
}

/* ************************************************************************* */
/**
 * @brief  Main control function for create optical depth diagnostics.
 *
 * @details
 *
 * This function is the main steering function for generating the optical depth
 * diagnostics.
 *
 * ************************************************************************** */

void
do_optical_depth_diagnostics (void)
{
  initialize_inclination_angles ();
  calculate_PI_edge_optical_depth ();
  create_optical_depth_spectrum ();
  free (INCLINATION_ANGLES);
}

/** ************************************************************************* */
/**
 * @file     py_optical_depth_output.c
 * @author   Edward Parkinson
 * @date     May 2021
 *
 * @brief    Functions for writing the output to stdout and to file.
 *
 * ************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "py_optd.h"

/* ************************************************************************* */
/**
 * @brief  Write a generic header to the provided file pointer.
 *
 * @details
 *
 * The version python, commit hash and date are written to the file pointer
 * provided. This assumes that the file pointer is pointing to the top of the
 * file, as this is supposed to be writing a header.
 *
 * ************************************************************************** */

void
write_generic_file_header (FILE *fp)
{
  char time_string[LINELENGTH];

  get_time (time_string);
  fprintf (fp, "# Python Version %s\n", VERSION);
  fprintf (fp, "# Git commit hash %s\n", GIT_COMMIT_HASH);
  fprintf (fp, "# Date %s\n", time_string);
  fprintf (fp, "#\n");
}

/* ************************************************************************** */
/**
 * @brief  Print the various optical depths for the photoionization edges
 *
 * @param[in]  SightLines_t  *inclinations  The inclination angles used
 * @param[in]  int  n_inclinations          The number of inclinations
 * @param[in]  Edges_t  edges[]             The photoionization edges
 * @param[in]  int  n_edges                 The number of edges evaluated
 * @param[in]  double *optical_depth        The optical depth of the edges
 * @param[in]  double *column_density       The column density of the inclinations
 *
 * @details
 *
 * Prints the different optical depths for each angle and optical depth.
 * Historically, this used to also print to file. But I found that it was mostly
 * useless to do this.
 *
 * ************************************************************************** */

int
print_optical_depths (SightLines_t *inclinations, int n_inclinations, PIEdge_t edges[], int n_edges, double *optical_depth,
                      double *column_density)
{
  printf ("Optical depths along the defined line of sights for domain %i:\n", CONFIG.domain);
  for (int i = 0; i < n_inclinations; i++)
  {
    printf ("%-8s:\n", inclinations[i].name);
    if (CONFIG.column_density == COLUMN_MODE_RHO)
    {
      printf (" -- Mass column density             : %3.2e cm^-2\n", column_density[i]);
      printf (" -- Approx. Hydrogen column density : %3.2e cm^-2\n", column_density[i] * rho2nh);
    }
    else
    {
      printf (" -- %s %i column density %-11s : %3.2e cm^-2\n", ele[ion[CONFIG.column_density_ion_number].nelem].name,
              ion[CONFIG.column_density_ion_number].istate, "", column_density[i]);
    }
    for (int j = 0; j < n_edges; j++)
    {
      printf (" -- %-31s : %3.2e\n", edges[j].name, optical_depth[i * n_edges + j]);
    }
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief  Write the optical depth spectrum to file.
 *
 * @param[in]  SightLines_t  *inclination  The inclinations angles the optical
 *                                         depth was extracted from.
 * @param[in]  int  n_inclinations         The number of inclination angles
 * @param[in]  double  *tau_spectrum       The optical depth spectrum values
 * @param[in]  double  freq_min            The starting frequency of the
 *                                         spectrum
 * @param[in]  double  dfreq               The frequency spacing of the spectrum
 *
 * @details
 *
 * Simply writes the optical depth spectra to the file named
 * root.tau_spec.
 *
 * ************************************************************************** */

int
save_optical_depth_spectrum (SightLines_t *inclinations, int n_inclinations, double *tau_spectrum, double freq_min, double d_freq)
{
  FILE *fp;
  char filename[LINELENGTH];

  int len = snprintf (filename, LINELENGTH, "%s.spec_tau", files.root);
  if (len < 0)
  {
    print_error ("error when creating filename string\n");
    return EXIT_FAILURE;
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    print_error ("unable to open %s in write mode\n", filename);
    return EXIT_FAILURE;
  }

  write_generic_file_header (fp);
  fprintf (fp, "%-15s %-15s ", "Freq.", "Lambda");
  for (int i = 0; i < n_inclinations; i++)
  {
    fprintf (fp, "%-15s ", inclinations[i].name);
  }
  fprintf (fp, "\n");

  double frequency = log10 (freq_min);
  for (int i = 0; i < NUM_FREQUENCY_BINS; i++)
  {
    double wavelength = VLIGHT / pow (10, frequency) / ANGSTROM;
    fprintf (fp, "%-15e %-15e ", pow (10, frequency), wavelength);

    for (int j = 0; j < n_inclinations; j++)
    {
      fprintf (fp, "%-15e ", tau_spectrum[j * NUM_FREQUENCY_BINS + i]);
    }

    fprintf (fp, "\n");
    frequency += d_freq;
  }

  if (fclose (fp))
  {
    print_error ("unable to close %s, output may be unfinished!\n", filename);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief  Write the photosphere location points to file.
 *
 * @param[in]  Positions_t *positions  An array of positions
 * @param[in]  int n_angles            The number of angles which were used
 *                                     in the photosphere calculation.
 *
 * @details
 *
 * This uses the Positions_t type which keeps track of the x, y and z locations
 * of where a photon reached a cumulative optical depth of TAU_DEPTH.
 *
 * ************************************************************************** */

void
write_photosphere_location_to_file (Pos_t *positions, int n_angles)
{
  int i;
  double pos1d[3];
  char filename[LINELENGTH];
  FILE *fp;

  int len = snprintf (filename, LINELENGTH, "%s.photosphere", files.root);
  if (len < 0)
  {
    print_error ("error when creating filename string\n");
    exit (EXIT_FAILURE);
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    print_error ("unable to open %s in write mode\n", filename);
    exit (EXIT_FAILURE);
  }

  write_generic_file_header (fp);
  fprintf (fp, "# Electron scatter photosphere locations for tau_es = %f\n#\n", CONFIG.arg_tau_surface);

  if (zdom[CONFIG.domain].coord_type != SPHERICAL)
  {
    fprintf (fp, "# %-15s %-15s %-15s $%-15s\n", "x", "y", "z", "rho");
  }
  else
  {
    fprintf (fp, "# %-15s\n", "r");
  }

  for (i = 0; i < n_angles; i++)
  {
    if (zdom[CONFIG.domain].coord_type != SPHERICAL)
    {
      double rho = sqrt (positions[i].x * positions[i].x + positions[i].y * positions[i].y);
      fprintf (fp, "%-15e %-15e %-15e %-15e\n", positions[i].x, positions[i].y, positions[i].z, rho);
    }
    else
    {
      pos1d[0] = positions[0].x;
      pos1d[1] = positions[0].y;
      pos1d[2] = positions[0].z;
      fprintf (fp, "%-15e\n", length (pos1d));
    }
  }

  if (fclose (fp))
  {
    print_error ("unable to close %s, output may be unfinished!\n", filename);
  }
}

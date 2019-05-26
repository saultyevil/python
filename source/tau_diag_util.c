/** ************************************************************************* */
/**
 * @file     tau_diag_util.c
 * @author   Edward Parkinson
 * @date     May 2019
 *
 * @brief    Contains utility functions for the optical depth integration
 *           algorithms - mostly for writing to file our stdout.
 *
 * @details
 *
 * ************************************************************************** */

#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "atomic.h"  // need otherwise it won't compile :^)
#include "python.h"  // for rank_global, LINELINE etc
#include "tau_diag.h"

/* ************************************************************************* */
/**
 * @brief           Print the various optical depths calculated using this
 *                  routine
 *
 * @param[in]       double   tau_store      The 2d array containing the optical
 *                                          depth for each observer and tau
 * @param[in]       double   col_den_store  The 2d array containing the column
 *                                          densities for each observer
 *
 * @return          void
 *
 * @details
 *
 * Simply prints the different optical depths for each angle and optical depth.
 * The same information will also be printed to an individual diag file named
 * optical_depth.diag which is located in the diag folder of the simulation.
 *
 * ************************************************************************** */

void
print_tau_table (double *tau_store, double *col_den_store)
{
  int itau;
  int ispec;
  int line_len;
  int const MAX_COL = 120;

  double tau;
  double col_den;

  char tmp_str[LINELENGTH];
  char observer_name[LINELENGTH];
  char diag_filename[3 * LINELENGTH];

  FILE *tau_diag;

  if (rank_global != 0)
    return;

  sprintf (diag_filename, "diag_%s/%s.optical_depth.diag", files.root, files.root);

  if (!(tau_diag = fopen (diag_filename, "w")))
  {
    Error ("%s:%s:%i: Unable to open optical depth diag file\n", __FILE__, __func__, __LINE__);
    Exit (1);
  }

  Log ("\nOptical depths along the defined line of sights:\n(-1 indicates an error occured)\n\n");
  fprintf (tau_diag, "Optical depths along the defined line of sights:\n(-1 indicates an error occured)\n\n");

  for (ispec = 0; ispec < N_ANGLES; ispec++)
  {
    strcpy (observer_name, tau_diag_observers[ispec].name);

    Log ("%s\n--------\n", observer_name);
    fprintf (tau_diag, "%s\n--------\n", observer_name);

    col_den = col_den_store[ispec];
    Log ("Column density: %3.2e g/cm^-2\n", col_den);
    fprintf (tau_diag, "Column density: %3.2e g/cm^-2\n", col_den);

    /*
     * For the photons which are pointing upwards and are not launched from the
     * origin
     */

    if (strcmp (observer_name, "A00P0.50") == 0)
    {
      Log ("Photon launched from x location: %e cm\n", print_xloc);
      fprintf (tau_diag, "Photon launched from x location: %e cm\n", print_xloc);
    }

    line_len = 0;
    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = tau_store[ispec * N_TAU + itau];
      line_len += sprintf (tmp_str, "tau_%-9s: %3.2e  ", TAU_DIAG_OPACS[itau].name, tau);

      if (line_len > MAX_COL)
      {
        line_len = 0;
        Log ("\n");
        fprintf (tau_diag, "\n");
      }

      Log ("%s", tmp_str);
      fprintf (tau_diag, "%s", tmp_str);
    }

    Log ("\n\n");
    fprintf (tau_diag, "\n\n");
  }

  if (fclose (tau_diag))
  {
    Error ("%s:%s:%i: could not close optical depth diag file\n", __FILE__, __func__, __LINE__);
    Exit (1);
  }
}

/* ************************************************************************* */
/**
 * @brief           Write the various optical depth spectra to file
 *
 * @param[in]       double tau_spectrum     The various optical depth spectra
 * @param[in]       double wave_min         The smallest wavelength in the
 *                                          spectra
 * @param[in]       double dwave            The size of the wavelength bins
 *
 * @return          void
 *
 * @details
 *
 * Simply write the optical depth spectra to the file named root.tau_spec.diag.
 * This file will be located in the diag folder.
 *
 * ************************************************************************** */

void
write_tau_spectrum (double *tau_spectrum, double wave_min, double dwave)
{
  int ispec;
  int iwave;

  double wavelength;

  char tau_spec_filename[3 * LINELENGTH];

  FILE *tau_spec_file;

  if (rank_global != 0)
    return;

  sprintf (tau_spec_filename, "diag_%s/%s.tau_spec.diag", files.root, files.root);

  if (!(tau_spec_file = fopen (tau_spec_filename, "w")))
  {
    Error ("%s:%s:%i: unable to open tau spectrum diag file\n", __FILE__, __func__, __LINE__);
    Exit (1);
  }

  /*
   * Write the file header
   */

  fprintf (tau_spec_file, "# Lambda ");
  for (ispec = 0; ispec < N_ANGLES; ispec++)
    fprintf (tau_spec_file, "%s ", tau_diag_observers[ispec].name);
  fprintf (tau_spec_file, "\n");

  /*
   * Write out the tau spectrum for each inclination angle
   */

  wavelength = wave_min;

  for (iwave = 0; iwave < NWAVE; iwave++)
  {
    fprintf (tau_spec_file, "%e ", wavelength);
    for (ispec = 0; ispec < N_ANGLES; ispec++)
      fprintf (tau_spec_file, "%e ", tau_spectrum[ispec * NWAVE + iwave]);
    fprintf (tau_spec_file, "\n");
    wavelength += dwave;
  }

  if (fclose (tau_spec_file))
  {
    Error ("%s:%s:%i: could not close tau spectrum diag file\n", __FILE__, __func__, __LINE__);
    Exit (1);
  }
}

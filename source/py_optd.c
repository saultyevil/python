/** ************************************************************************* */
/**
 * @file  py_optical_depth.c
 * @author  Edward Parkinson
 * @date  February 2021
 *
 * @brief  File containing the main functions defining the program.
 *
 * ************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include "atomic.h"
#include "python.h"
#include "py_optd.h"

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

void
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
    print_error ("Unable to allocate memory for domain\n");
    exit (EXIT_FAILURE);
  }

  if (wind_read (windsave_filename) < 0)
  {
    print_error ("unable to open %s\n", windsave_filename);
    exit (EXIT_FAILURE);
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
        print_error ("Unable to open %s when spectrum cycles have been run\n", specsave_filename);
        exit (EXIT_FAILURE);
      }
    }
  }
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
calculate_photoionization_optical_depths (double *input_inclinations)
{
  int i, j;
  int err;
  double c_frequency, c_optical_depth, c_column_density;
  double *optical_depth_values = NULL, *column_density_values = NULL;
  struct photon photon;
  enum RunModeEnum original_run_mode = RUN_MODE;
  RUN_MODE = MODE_IGNORE_ELECTRON_SCATTERING;

  Edges_t edges[] = {
    {"HLymanEdge", 3.387485e+15},
    {"HBalmerEdge", 8.293014e+14},
    {"HeI24eVEdge", 5.9483e+15},
    {"HeII54eVEdge", 1.394384e+16}
  };

  const int n_edges = sizeof edges / sizeof edges[0];

  int n_inclinations;
  SightLines_t *inclinations = initialize_inclination_angles (&n_inclinations, input_inclinations);

  optical_depth_values = calloc (n_inclinations * n_edges, sizeof *optical_depth_values);
  if (optical_depth_values == NULL)
  {
    print_error ("cannot allocate %lu bytes for optical_depths\n", n_inclinations * n_edges * sizeof *optical_depth_values);
    exit (EXIT_FAILURE);
  }

  column_density_values = calloc (n_inclinations, sizeof *column_density_values);
  if (column_density_values == NULL)
  {
    print_error ("cannot allocate %lu bytes for column_densities\n", n_inclinations * sizeof *column_density_values);
    exit (EXIT_FAILURE);
  }

  /*
   * Now extract the optical depths and mass column densities. We loop over
   * each PI edge for each inclination angle.
   */

  for (i = 0; i < n_inclinations; i++)
  {
    for (j = 0; j < n_edges; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;
      c_frequency = edges[j].freq;

      err = create_photon (&photon, c_frequency, inclinations[i].lmn);
      if (err == EXIT_FAILURE)
      {
        print_error ("skipping photon of frequency %e\n", c_frequency);
        continue;
      }

      err = integrate_tau_across_wind (&photon, &c_column_density, &c_optical_depth);
      if (err == EXIT_FAILURE)
        continue;               // do not throw extra warning when one is already thrown in integrate_tau_across_wind

      optical_depth_values[i * n_edges + j] = c_optical_depth;
      column_density_values[i] = c_column_density;
    }
  }

  print_optical_depths (inclinations, n_inclinations, edges, n_edges, optical_depth_values, column_density_values);
  free (inclinations);
  free (optical_depth_values);
  free (column_density_values);
  RUN_MODE = original_run_mode;
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
 * from this location towards the observer where it escapes, where integrate_tau_across_wind
 * returns the integrated optical depth along its path to escape. This is done
 * for a range of photon frequencies to find how optical depth changes with
 * frequency.
 *
 * This processes can take some time compared to tau_evalulate_photo_edges. But,
 * since NUM_FREQUENCY_BINS photons are being generated for each spectrum and the fact
 * that these photons do not interact, the spectra does not actually take that
 * long to complete.
 *
 * ************************************************************************** */

void
create_optical_depth_spectrum (double u_freq_min, double u_freq_max, double *input_inclinations)
{
  int i, j;
  int err;
  double *tau_spectrum;
  double c_optical_depth, c_column_density;
  double c_frequency, freq_min, freq_max, d_freq;
  struct photon photon;

  int n_inclinations;
  SightLines_t *inclinations = initialize_inclination_angles (&n_inclinations, input_inclinations);

  printf ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (n_inclinations * NUM_FREQUENCY_BINS, sizeof *tau_spectrum);
  if (tau_spectrum == NULL)
  {
    print_error ("cannot allocate %lu bytes for tau_spectrum\n", n_inclinations * NUM_FREQUENCY_BINS * sizeof *tau_spectrum);
    exit (EXIT_FAILURE);
  }

  /*
   * We have a complicated if statement first, though. If a freq_min
   * or a freq_max was provided, then we need to get this first and set
   * the frequency limits appropriately. If neither are defined, then we will
   * use some hardwired limits. The frequency range of the extracted will be
   * used, however if xxpsec is NULL (no observer spectrum exists), then the
   * frequency range will be over a default 100 - 10,000 Angstrom band.
   */

  if (u_freq_min > 0 || u_freq_max > 0)
  {
    if (u_freq_min > 0)
    {
      freq_min = u_freq_min;
    }
    else
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
    }

    if (u_freq_max > 0)
    {
      freq_max = u_freq_max;
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

  d_freq = (log10 (freq_max) - log10 (freq_min)) / NUM_FREQUENCY_BINS;
  kbf_need (freq_min, freq_max);

  /*
   * Now create the optical depth spectra for each inclination
   */

  for (i = 0; i < n_inclinations; i++)
  {
    printf ("  - Creating spectrum: %s\n", inclinations[i].name);
    c_frequency = log10 (freq_min);

    for (j = 0; j < NUM_FREQUENCY_BINS; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;

      err = create_photon (&photon, pow (10, c_frequency), inclinations[i].lmn);
      if (err == EXIT_FAILURE)
      {
        print_error ("skipping photon of frequency %e\n", pow (10, c_frequency));
        continue;
      }

      err = integrate_tau_across_wind (&photon, &c_column_density, &c_optical_depth);
      if (err == EXIT_FAILURE)
        continue;

      tau_spectrum[i * NUM_FREQUENCY_BINS + j] = c_optical_depth;
      c_frequency += d_freq;
    }
  }

  write_optical_depth_spectrum (inclinations, n_inclinations, tau_spectrum, freq_min, d_freq);
  free (tau_spectrum);
  free (inclinations);
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

void
calculate_cell_optical_depth_spectrum (double u_freq_min, double u_freq_max, double *input_inclinations)
{
  int i, j;
  int err;
  double *tau_spectrum;
  double c_optical_depth, c_column_density;
  double c_frequency, freq_min, freq_max, d_freq;
  struct photon photon;

  int n_inclinations;
  SightLines_t *inclinations = initialize_inclination_angles (&n_inclinations, input_inclinations);

  printf ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (n_inclinations * NUM_FREQUENCY_BINS, sizeof *tau_spectrum);
  if (tau_spectrum == NULL)
  {
    print_error ("cannot allocate %lu bytes for tau_spectrum\n", n_inclinations * NUM_FREQUENCY_BINS * sizeof *tau_spectrum);
    exit (EXIT_FAILURE);
  }

  /*
   * We have a complicated if statement first, though. If a freq_min
   * or a freq_max was provided, then we need to get this first and set
   * the frequency limits appropriately. If neither are defined, then we will
   * use some hardwired limits. The frequency range of the extracted will be
   * used, however if xxpsec is NULL (no observer spectrum exists), then the
   * frequency range will be over a default 100 - 10,000 Angstrom band.
   */

  if (u_freq_min > 0 || u_freq_max > 0)
  {
    if (u_freq_min > 0)
    {
      freq_min = u_freq_min;
    }
    else
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
    }

    if (u_freq_max > 0)
    {
      freq_max = u_freq_max;
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

  d_freq = (log10 (freq_max) - log10 (freq_min)) / NUM_FREQUENCY_BINS;
  kbf_need (freq_min, freq_max);

  /*
   * Now create the optical depth spectra for each inclination
   */

  for (i = 0; i < n_inclinations; i++)
  {
    printf ("  - Creating spectrum: %s\n", inclinations[i].name);
    c_frequency = log10 (freq_min);

    for (j = 0; j < NUM_FREQUENCY_BINS; j++)
    {
      c_optical_depth = 0.0;

      err = create_photon (&photon, pow (10, c_frequency), inclinations[i].lmn);
      if (err == EXIT_FAILURE)
      {
        print_error ("skipping photon of frequency %e\n", pow (10, c_frequency));
        continue;
      }

      err = optical_depth_across_cell (&photon, &c_optical_depth);
      if (err == EXIT_FAILURE)
        continue;

      tau_spectrum[i * NUM_FREQUENCY_BINS + j] = c_optical_depth;
      c_frequency += d_freq;
    }
  }

  printf ("huh???\n");
  write_optical_depth_spectrum (inclinations, n_inclinations, tau_spectrum, freq_min, d_freq);
  free (tau_spectrum);
  free (inclinations);
}

/* ************************************************************************* */
/**
 * @brief   Find the electron scattering photosphere.
 *
 * @details
 *
 * This is the main controlling function for finding the photosphere. The
 * electron scattering optical depth is controlled by the global variable
 * TAU_DEPTH.
 *
 * ************************************************************************** */

void
calculate_optical_depth_surfaces (void)
{
  int i, error;
  double optical_depth, column_density;
  struct photon photon;
  SightLines_t *inclinations;
  double input_inclinations[MAX_CUSTOM_ANGLES];

  for (i = 0; i < MAX_CUSTOM_ANGLES; ++i)       // required for initialize_inclination_angles, even though unused
    input_inclinations[i] = -1.0;

  int n_inclinations;
  inclinations = initialize_inclination_angles (&n_inclinations, input_inclinations);

  Positions_t *positions = calloc (n_inclinations, sizeof (Positions_t));
  if (positions == NULL)
  {
    print_error ("unable to allocate memory for the positions array\n");
    exit (EXIT_FAILURE);
  }

  const double test_freq = 8e14;        // todo: this probably need to be a possible input

  printf ("Locating electron scattering photosphere surface for tau_es = %f\n", TAU_DEPTH);

  for (i = 0; i < n_inclinations; i++)
  {
    error = create_photon (&photon, test_freq, inclinations[i].lmn);
    if (error)
    {
      positions[i].x = positions[i].y = positions[i].z = -1.0;
      continue;
    }

    optical_depth = column_density = 0;

    error = integrate_tau_across_wind (&photon, &column_density, &optical_depth);
    if (error)
    {
      positions[i].x = positions[i].y = positions[i].z = -1.0;
      continue;
    }

    positions[i].x = photon.x[0];
    positions[i].y = photon.x[1];
    positions[i].z = photon.x[2];
  }

  write_photosphere_location_to_file (positions, n_inclinations);
  free (inclinations);
  free (positions);
}

/* ************************************************************************* */
/**
 * @brief  The main function of the program.
 *
 * @param  argc  The number of command line arguments
 * @param  argv  The command line arguments
 *
 * @return  EXIT_SUCCESS
 *
 * @details
 *
 * ************************************************************************** */

int
main (int argc, char *argv[])
{
  timer ();
  Log_set_verbosity (2);
  Log_print_max (10);
  Log_quit_after_n_errors ((int) 1e8);
  init_rand ((int) time (NULL));

  rel_mode = REL_MODE_FULL;
  SMAX_FRAC = 1.0;
  DENSITY_PHOT_MIN = 1.e-10;
  COLUMN_MODE = COLUMN_MODE_RHO;
  RUN_MODE = MODE_SPECTRUM;
  N_DOMAIN = 0;

  struct CommandlineArguments arguments;
  arguments = get_arguments (argc, argv);


  printf ("%-20s Optical depth diagnostics beginning\n", "TAU");
  initialize_wind_structures ();

  switch (RUN_MODE)
  {
  case MODE_PHOTOION:
    calculate_photoionization_optical_depths (arguments.inclinations);
    break;
  case MODE_SPECTRUM:
    create_optical_depth_spectrum (arguments.freq_min, arguments.freq_max, arguments.inclinations);
    break;
  case MODE_CELL_SPECTRUM:
    calculate_cell_optical_depth_spectrum (arguments.freq_min, arguments.freq_max, arguments.inclinations);
    break;
  case MODE_SURFACE:
    SMAX_FRAC = 0.1;            // takes longer, but higher accuracy
    calculate_optical_depth_surfaces ();
    break;
  default:
    print_error ("Mode %d is an unknown run mode, not sure how you got here so exiting the program\n", RUN_MODE);
    break;
  }

  printf ("\n%-20s Optical depth diagnostics completed\n", "TAU");
  printf ("Completed optical depth diagnostics. The elapsed TIME was %f\n", timer ());
  error_summary ("end of program");

  return EXIT_SUCCESS;
}

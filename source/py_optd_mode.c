/** ************************************************************************* */
/**
 * @file  py_optd_mode.c
 * @author  Edward Parkinson
 * @date  September 2024
 *
 * @brief
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "atomic.h"
#include "python.h"
#include "py_optd.h"

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

int
calculate_photoionization_optical_depths (void)
{
  int error;
  struct photon photon;
  struct SightLines inclinations[MAX_ANGLES];   // cba with mallocing this
  CONFIG.ignore_electron_scattering = true;

  PIEdge_t pi_edges[] = {
    {"H_Lyman", 3.387485e+15},
    {"H_Balmer", 8.293014e+14},
    {"HeI_24eV", 5.948300e+15},
    {"HeII_54eV", 1.394384e+16}
  };
  const int n_edges = sizeof pi_edges / sizeof pi_edges[0];

  if (CONFIG.column_density_ion_number > nions - 1)
  {
    printf ("The ion number %i is an invalid ion number: there are %i ions available\n", CONFIG.column_density_ion_number, nions);
    return EXIT_FAILURE;
  }
  if (CONFIG.column_density == COLUMN_MODE_ION)
  {
    printf ("Extracting column density for %s %i\n", ele[ion[CONFIG.column_density_ion_number].nelem].name,
            ion[CONFIG.column_density_ion_number].istate);
  }

  initialize_inclination_angles (inclinations);
  double *optical_depths = calloc (CONFIG.arg_num_inc * n_edges, sizeof *optical_depths);
  if (optical_depths == NULL)
  {
    print_error ("failed to allocate %lu bytes for optical depths\n", CONFIG.arg_num_inc * n_edges * sizeof *optical_depths);
    return EXIT_FAILURE;
  }
  double *column_densities = calloc (CONFIG.arg_num_inc, sizeof *column_densities);
  if (column_densities == NULL)
  {
    print_error ("failed to allocate %lu bytes for column densities\n", CONFIG.arg_num_inc * sizeof *column_densities);
    return EXIT_FAILURE;
  }

  /*
   * Now extract the optical depths and mass column densities. We loop over
   * each PI edge for each inclination angle.
   */

  for (int i = 0; i < CONFIG.arg_num_inc; i++)
  {
    for (int j = 0; j < n_edges; j++)
    {
      double optical_depth = 0.0;
      double column_density = 0.0;
      const double frequency = pi_edges[j].freq;

      initialise_photon_packet (&photon, frequency, inclinations[i].direction_vector);

      error = integrate_tau_across_wind (&photon, &column_density, &optical_depth);
      if (error)
      {
        continue;               // we can deal with not finding an optical depth here
      }

      optical_depths[i * n_edges + j] = optical_depth;
      column_densities[i] = column_density;
    }
  }

  error = print_optical_depths (inclinations, CONFIG.arg_num_inc, pi_edges, n_edges, optical_depths, column_densities);

  free (optical_depths);
  free (column_densities);

  return error;
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

int
create_optical_depth_spectrum (void)
{
  int error;
  struct photon photon;

  SightLines_t inclinations[MAX_ANGLES];

  initialize_inclination_angles (inclinations);
  int n_inclinations = CONFIG.arg_num_inc;      // TODO: this is nasty!
  double *spectrum = calloc (n_inclinations * NUM_FREQUENCY_BINS, sizeof (double));
  if (spectrum == NULL)
  {
    print_error ("Unable to allocate %lu bytes for spectrum\n", n_inclinations * NUM_FREQUENCY_BINS * sizeof (double));
    exit (EXIT_FAILURE);
  }

  // Deal with frequency range
  set_frequency_range ();       // TODO: don't like doing this with global state!!
  const double freq_max = CONFIG.arg_freq_max;
  const double freq_min = CONFIG.arg_freq_min;
  const double d_freq = (log10 (freq_max) - log10 (freq_min)) / NUM_FREQUENCY_BINS;
  const int output_freq = NUM_FREQUENCY_BINS / 10;

  // This is an optimisation to reduce the number of opacities we need to check
  kbf_need (freq_min, freq_max);

  // Calculate optical depth spectrum for each inclination
  printf ("Creating optical depth spectra between %e Hz - %e Hz for %d observer angles\n", freq_min, freq_max, n_inclinations);
  for (int i = 0; i < n_inclinations; i++)
  {
    printf (" -- %s ", inclinations[i].name);
    double frequency = log10 (freq_min);
    for (int j = 0; j < NUM_FREQUENCY_BINS; j++)
    {
      if (j % output_freq == 0)
      {
        printf (".");
        fflush (stdout);
      }

      double optical_depth = 0.0;
      double column_density = 0.0;

      error = initialise_photon_packet (&photon, pow (10, frequency), inclinations[i].direction_vector);
      if (error == EXIT_FAILURE)
      {
        print_error ("Issue for photon for frequency bin %e Hz\n", pow (10, frequency));
        continue;
      }

      error = integrate_tau_across_wind (&photon, &column_density, &optical_depth);
      if (error == EXIT_FAILURE)
      {
        continue;
      }

      spectrum[i * NUM_FREQUENCY_BINS + j] = optical_depth;
      frequency += d_freq;
    }
    printf ("\n");
  }

  error = save_optical_depth_spectrum (inclinations, n_inclinations, spectrum, freq_min, d_freq);
  if (error)
  {
    print_error ("Unable to write optical depth spectrum to disk\n");
  }
  free (spectrum);

  return error ? EXIT_FAILURE : EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

int
create_cell_optical_depth_spectrum (void)
{
  // int i, j;
  // int err;
  // double *tau_spectrum;
  // double c_optical_depth, c_column_density;
  // double c_frequency, freq_min, freq_max, d_freq;
  // struct photon photon;

  // SightLines_t inclinations[MAX_ANGLES];

  // initialize_inclination_angles(inclinations);
  // int n_inclinations = CONFIGURATION.n_inclinations;

  // printf("Creating optical depth spectra:\n");

  // tau_spectrum = calloc(n_inclinations * NUM_FREQUENCY_BINS, sizeof *tau_spectrum);
  // if (tau_spectrum == NULL)
  // {
  //     print_error("cannot allocate %lu bytes for tau_spectrum\n", n_inclinations * NUM_FREQUENCY_BINS * sizeof *tau_spectrum);
  //     exit(EXIT_FAILURE);
  // }

  // /*
  //  * We have a complicated if statement first, though. If a freq_min
  //  * or a freq_max was provided, then we need to get this first and set
  //  * the frequency limits appropriately. If neither are defined, then we will
  //  * use some hardwired limits. The frequency range of the extracted will be
  //  * used, however if xxpsec is NULL (no observer spectrum exists), then the
  //  * frequency range will be over a default 100 - 10,000 Angstrom band.
  //  */

  // if (u_freq_min > 0 || u_freq_max > 0)
  // {
  //     if (u_freq_min > 0)
  //     {
  //         freq_min = u_freq_min;
  //     }
  //     else
  //     {
  //         freq_min = VLIGHT / (10000 * ANGSTROM);
  //     }

  //     if (u_freq_max > 0)
  //     {
  //         freq_max = u_freq_max;
  //     }
  //     else
  //     {
  //         freq_max = VLIGHT / (100 * ANGSTROM);
  //     }

  //     if (freq_max < freq_min)
  //     {
  //         print_error("frequency range given has set freq_max (%e) < freq_min (%e) \n", freq_max, freq_min);
  //         exit(EXIT_FAILURE);
  //     }
  // }
  // else
  // {
  //     if ((geo.nangles == 0 && xxspec == NULL) || (geo.swavemax == 0 && geo.swavemin == 0))
  //     {
  //         freq_min = VLIGHT / (10000 * ANGSTROM);
  //         freq_max = VLIGHT / (100 * ANGSTROM);
  //     }
  //     else
  //     {
  //         freq_min = VLIGHT / (geo.swavemax * ANGSTROM);
  //         freq_max = VLIGHT / (geo.swavemin * ANGSTROM);
  //         if (sane_check(freq_min))
  //         {
  //             freq_min = VLIGHT / (10000 * ANGSTROM);
  //             print_error("freq_min has an invalid value setting to %e\n", freq_min);
  //         }
  //         if (sane_check(freq_max))
  //         {
  //             freq_max = VLIGHT / (100 * ANGSTROM);
  //             print_error("freq_min has an invalid value setting to %e\n", freq_max);
  //         }
  //     }
  // }

  // d_freq = (log10(freq_max) - log10(freq_min)) / NUM_FREQUENCY_BINS;
  // kbf_need(freq_min, freq_max);

  // /*
  //  * Now create the optical depth spectra for each inclination
  //  */

  // for (i = 0; i < n_inclinations; i++)
  // {
  //     printf("  - Creating spectrum: %s\n", inclinations[i].name);
  //     c_frequency = log10(freq_min);

  //     for (j = 0; j < NUM_FREQUENCY_BINS; j++)
  //     {
  //         c_optical_depth = 0.0;

  //         err = initialise_photon_packet(&photon, pow(10, c_frequency), inclinations[i].direction_vector);
  //         if (err == EXIT_FAILURE)
  //         {
  //             print_error("skipping photon of frequency %e\n", pow(10, c_frequency));
  //             continue;
  //         }

  //         err = move_photon_across_cell(&photon, &c_optical_depth);
  //         if (err == EXIT_FAILURE)
  //             continue;

  //         tau_spectrum[i * NUM_FREQUENCY_BINS + j] = c_optical_depth;
  //         c_frequency += d_freq;
  //     }
  // }

  // write_optical_depth_spectrum(inclinations, n_inclinations, tau_spectrum, freq_min, d_freq);
  // free(tau_spectrum);
  // free(inclinations);

  return EXIT_SUCCESS;
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

int
find_optical_depth_surface (void)
{
  struct photon photon;
  struct SightLines inclinations[MAX_ANGLES];
  const double TEST_FREQUENCY = 8e14;   // TODO: this probably need to be a possible input

  SMAX_FRAC = 0.01;             // takes longer, but higher accuracy
  initialize_inclination_angles (inclinations);

  Pos_t *positions = calloc (CONFIG.arg_num_inc, sizeof (Pos_t));
  if (positions == NULL)
  {
    print_error ("Unable to allocate memory to store coordinates of optical depth surface\n");
    return EXIT_FAILURE;
  }

  printf ("Locating electron scattering surface of constant optical depth (to escape) for tau_es = %f\n", CONFIG.arg_tau_surface);
  const int OUTPUT_FREQ = CONFIG.arg_num_inc / 10;

  for (int i = 0; i < CONFIG.arg_num_inc; i++)
  {
    if (i % OUTPUT_FREQ == 0)
    {
      printf (".");
      fflush (stdout);
    }

    initialise_photon_packet (&photon, TEST_FREQUENCY, inclinations[i].direction_vector);
    double optical_depth = 0;
    double column_density = 0;
    const int error = integrate_tau_across_wind (&photon, &column_density, &optical_depth);
    if (error)
    {
      print_error ("Error occurred -- returning without finishing");
      return EXIT_FAILURE;
    }

    positions[i].x = photon.x[0];
    positions[i].y = photon.x[1];
    positions[i].z = photon.x[2];
  }

  printf ("\n");
  write_photosphere_location_to_file (positions, CONFIG.arg_num_inc);
  free (positions);

  return EXIT_SUCCESS;
}

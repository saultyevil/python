/** ************************************************************************* */
/**
 * @file py_optd.c
 * @author Edward Parkinson
 * @date February 2021
 *
 * @brief File containing the main functions defining the program.
 *
 * ************************************************************************** */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "atomic.h"
#include "python.h"
#include "py_optd.h"

/* ************************************************************************* */
/**
 * @brief Set default values for the configuration of the program.
 *
 * @details
 *
 * Populates values for the CONFIG struct and initialises some global values
 * used in Sirocco/Python during photon transport and logging.
 *
 * ************************************************************************** */

void
set_default_configuration (void)
{
  // Initialise things to make Sirocco/Python work
  // In the future, some may need to be changed via CLI arguments
  Log_set_verbosity (1);
  Log_print_max (10);
  Log_quit_after_n_errors ((int) 1e8);
  init_rand ((int) time (NULL));
  rel_mode = REL_MODE_FULL;
  SMAX_FRAC = 1.0;
  DENSITY_PHOT_MIN = 1.e-10;

  // Initialise basic controlling parameters
  CONFIG.mode = MODE_SPECTRUM;
  CONFIG.column_density = COLUMN_MODE_RHO;
  CONFIG.column_density_ion_number = 0;
  CONFIG.domain = 0;
  CONFIG.ignore_electron_scattering = false;

  // Initialise values for CLI arguments
  CONFIG.arg_freq_min = -1.0;
  CONFIG.arg_freq_max = -1.0;
  CONFIG.arg_num_inc = 0;
  CONFIG.arg_tau_surface = 0.0;
  for (int i = 0; i < MAX_ANGLES; i++)
  {
    CONFIG.arg_inclinations[i] = -1.0;
  }
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
  int error = EXIT_SUCCESS;

  timer ();
  set_default_configuration ();
  parse_optd_arguments (argc, argv);
  printf ("%-20s Optical depth diagnostics beginning\n", "TAU");
  error = initialize_wind_structures ();

  // Bit of an dirty way to deal with errors when initialising the simulation
  // data -- but it will keep the output consistent.
  if (error != EXIT_FAILURE)
  {
    switch (CONFIG.mode)
    {
    case MODE_PHOTOION:
      error = calculate_photoionization_optical_depths ();
      break;
    case MODE_SPECTRUM:
      error = create_optical_depth_spectrum ();
      break;
    case MODE_CELL_SPECTRUM:
      error = create_cell_optical_depth_spectrum ();
      break;
    case MODE_SURFACE:
      error = find_optical_depth_surface ();
      break;
    }
  }

  if (error)
  {
    printf ("%-20s Optical depth diagnostics failed\n", "TAU");
    printf ("Optical depth diagnostics failed. The elapsed TIME was %f\n", timer ());
  }
  else
  {
    printf ("%-20s Optical depth diagnostics completed\n", "TAU");
    printf ("Completed optical depth diagnostics. The elapsed TIME was %f\n", timer ());
  }

  return error;
}

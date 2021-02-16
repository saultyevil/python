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
#include <unistd.h>

#include "atomic.h"
#include "python.h"
#include "py_optical_depth.h"

/* ************************************************************************* */
/**
 * @brief  Help text for the program.
 *
 * ************************************************************************** */

void
print_help (void)
{
  char *help = "Create optical depth diagnostics for a Python simulation.\n\n"
    "usage: py_optical_depth [-h] [-cion nion] [-classic] [--version] root\n"
    "\n"
    "-h           Print this help message and exit\n"
    "-cion nion   Extract the column density for an ion of number nion\n"
    "-classic     Use linear frequency transforms. Use when Python run in classic mode.\n"
    "--version    Print the version information and exit.\n";
  printf ("%s", help);
}

/* ************************************************************************* */
/**
 * @brief  Parse run time arguments provided at the command line.
 *
 * @param[in]  argc  The number of arguments provided
 * @param[in]  argv  The command line arguments
 *
 * @details
 *
 * Reads the command line arguments. Assumes the final argument is the root name
 * of the model. If no arguments are provided, then the program will run in a
 * default mode and query the user for the root name of the simulation.
 *
 * ************************************************************************** */

void
get_arguments (int argc, char *argv[])
{
  int n_read;
  char input[LINELENGTH];

  /*
   * If no command line arguments are provided, then use fgets to query the
   * root name from the user and return
   */

  if (argc == 1)
  {
    printf ("Please input the model root name: ");
    if (fgets (input, LINELENGTH, stdin) == NULL)
    {
      printf ("Unable to parse root name\n");
      exit (EXIT_FAILURE);
    }
    get_root (files.root, input);
    return;
  }

  /*
   * Otherwise, iterate over the command line arguments and initialize various
   * variables
   */

  n_read = 0;
  for (int i = 1; i < argc; ++i)
  {
    if (!strcmp (argv[i], "-h"))
    {
      print_help ();
      exit (EXIT_SUCCESS);
    }
    else if (!strcmp (argv[i], "-classic"))
    {
      rel_mode = REL_MODE_LINEAR;
      n_read = i;
    }
    else if (!strcmp (argv[i], "--version"))
    {
      printf ("Python version %s\n", VERSION);
      printf ("Built from git commit %s\n", GIT_COMMIT_HASH);
      if (GIT_DIFF_STATUS)
        printf ("This version was compiled with %d files with uncommited changes\n", GIT_DIFF_STATUS);
      exit (EXIT_SUCCESS);
    }
    else if (!strcmp (argv[i], "-cion"))
    {
      char *check;
      COLUMN_MODE = COLUMN_MODE_ION;
      COLUMN_MODE_ION_NUMBER = (int) strtol (argv[i + 1], &check, 10);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -cion into an integer\n");
        exit (EXIT_FAILURE);
      }
      if (COLUMN_MODE_ION_NUMBER < 0)
      {
        printf ("Argument for -cion cannot be negative\n");
        exit (EXIT_FAILURE);
      }
      i++;
      n_read = i;
    }
    else if (!strncmp (argv[i], "-", 1))
    {
      printf ("Unknown argument %s\n", argv[i]);
      print_help ();
      exit (EXIT_FAILURE);
    }
  }

  if (n_read + 1 == argc)
  {
    printf ("All command line arguments have been consumed without specifying a root name\n");
    exit (EXIT_FAILURE);
  }

  get_root (files.root, argv[argc - 1]);
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
  char windsave_filename[LINELENGTH + 24];
  char specsave_filename[LINELENGTH + 24];
  void do_optical_depth_diagnostics (void);

  timer ();

  /*
   * Initialize global variables which are required for photon
   * transport and general program operation
   */

  Log_set_verbosity (0);
  init_rand (time (NULL));
  rel_mode = REL_MODE_FULL;     // this is updated in get_arguments if required
  SMAX_FRAC = 0.5;
  DENSITY_PHOT_MIN = 1.e-10;
  COLUMN_MODE = COLUMN_MODE_RHO;
  N_INCLINATION_ANGLES = 0;

  get_arguments (argc, argv);
  printf ("%-20s Optical depth diagnostics beginning\n", "TAU");
  strcpy (windsave_filename, files.root);
  strcat (windsave_filename, ".wind_save");
  strcpy (specsave_filename, files.root);
  strcat (specsave_filename, ".spec_save");

  /*
   * Read in the wind_save file and initialize the wind cones and DFUDGE which
   * are important for photon transport. The atomic data is also read in at
   * this point (which is also very important)
   */

  if (wind_read (windsave_filename) < 0)
  {
    printf ("Unable to open %s\n", windsave_filename);
    exit (EXIT_FAILURE);
  }

  DFUDGE = setup_dfudge ();
  setup_windcone ();

  /*
   * Do some further error checking on the COLUMN_MODE_ION to ensure that
   * it is a sensible number and print the ion to be extracted
   */

  if (COLUMN_MODE == COLUMN_MODE_ION)
  {
    if (COLUMN_MODE_ION_NUMBER > nions - 1)
    {
      printf ("Argument for -cion > nions\n");
      exit (1);
    }
    printf ("Extracting column density for %s %i\n", ele[ion[COLUMN_MODE_ION_NUMBER].nelem].name, ion[COLUMN_MODE_ION_NUMBER].istate);
  }
  else
  {
    printf ("Extracting mass and Hydrogen column density\n");
  }

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
        printf ("Unable to open %s when spectrum cycles have been run\n", specsave_filename);
        exit (EXIT_FAILURE);
      }
    }
  }

  /*
   * Now we can finally being the optical depth diagnostics...
   */

  do_optical_depth_diagnostics ();
  printf ("\n%-20s Optical depth diagnostics completed\n", "TAU");
  printf ("Completed optical depth diagnostics. The elapsed TIME was %f\n", timer ());

  return EXIT_SUCCESS;
}

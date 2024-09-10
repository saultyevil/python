/** ************************************************************************* */
/**
 * @file  py_optd_parse.c
 * @author Edward Parkinson
 * @date August 2024
 *
 * ************************************************************************** */

#include <stdbool.h>
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
 * @brief  Help text for the program.
 *
 * ************************************************************************** */

static void
print_help (void)
{
  char *help = "A utility program to analyse the optical depth in a Python model.\n\n"
    "usage: py_otpd pi|wind|cell|surface [arguments] root\n";
  printf ("%s", help);
}

/* ************************************************************************* */
/**
 * @brief Parse arguments for mode 'pi'
 *
 * @param[in]  argc  The number of arguments provided
 * @param[in]  argv  The command line arguments
 *
 * @return num_args_processed  The number of arguments processed
 *
 * @details
 *
 * ************************************************************************** */

static int
parse_args_for_mode_pi (int argc, char *argv[])
{
  char *ret;
  int num_args_processed = 0;

  // If --ion is used, then we extract an ion's column density rather than
  // the hydrogen column density
  if (strcmp (argv[2], "-ion") == 0)
  {
    if (argc < 4)
    {
      printf ("Invalid number of arguments for mode 'pi -ion'\n");
      print_help ();
      exit (EXIT_FAILURE);
    }
    CONFIG.column_density = COLUMN_MODE_ION;
    CONFIG.column_density_ion_number = (int) strtol (argv[3], &ret, 10);
    if (*ret != '\0')
    {
      printf ("Unable to convert argument provided for -ion into an integer\n");
      exit (EXIT_FAILURE);
    }
    if (CONFIG.column_density_ion_number < 0)
    {
      printf ("Invalid choice for -ion!\n");
      print_help ();
      exit (EXIT_FAILURE);
    }

    num_args_processed = 2;
  }

  return num_args_processed;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]  argc  The number of arguments provided
 * @param[in]  argv  The command line arguments
 *
 * @details
 *
 * ************************************************************************** */

static int
parse_args_for_mode_wind (int argc, char *argv[])
{
  int num_args_processed = 0;

  for (int i = 2; i < argc; ++i)
  {
    if (!strcmp (argv[i], "-i"))
    {
      char *end_ptr;

      CONFIG.arg_num_inc = 0;
      for (int j = 1; j < MAX_CUSTOM_ANGLES + 1; ++j)
      {
        float inclination = (float) strtod (argv[i + j], &end_ptr);
        if (inclination < 0 || inclination > 90)
        {
          printf ("Invalid inclination %f°! Has to be 0° <= inclination <= 90°\n", inclination);
          exit (EXIT_FAILURE);
        }
        if (strcmp (argv[i + j], end_ptr) == 0) // conversion error indicates no more numbers
        {
          break;
        }
        CONFIG.arg_inclinations[j - 1] = inclination;
        CONFIG.arg_num_inc += 1;
      }

      if (CONFIG.arg_num_inc == 0)
      {
        printf ("Invalid input: no inclination were provided for --i option\n");
        exit (EXIT_FAILURE);
      }

      num_args_processed += 1 + CONFIG.arg_num_inc;
      i += num_args_processed;  // avoid trying to reprocess the inclination angles
    }
  }

  return num_args_processed;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]  argc  The number of arguments provided
 * @param[in]  argv  The command line arguments
 *
 * @details
 *
 * ************************************************************************** */

static int
parse_args_for_mode_cell (int argc, char *argv[])
{
  char *ret;
  int num_args_processed = 0;

  if (argc < 4)
  {
    printf ("Invalid arguments for mode 'cell'!\n");
    print_help ();
    exit (EXIT_FAILURE);
  }

  // Get i and j coordinates for cell in question
  int i = (int) strtol (argv[2], &ret, 10);
  if (*ret != '\0')
  {
    printf ("Unable to convert argument for i-th coordinate into an integer\n");
    exit (EXIT_FAILURE);
  }
  int j = (int) strtol (argv[3], &ret, 10);
  if (*ret != '\0')
  {
    printf ("Unable to convert argument for j-th coordinate into an integer\n");
    exit (EXIT_FAILURE);
  }
  num_args_processed += 2;

  // We can't do anything with this yet, as we don't know the number of elements
  // in the wind
  CONFIG.arg_wind_i = i;
  CONFIG.arg_wind_j = j;

  return num_args_processed;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]  argc  The number of arguments provided
 * @param[in]  argv  The command line arguments
 *
 * @details
 *
 * ************************************************************************** */

static int
parse_args_for_mode_surface (int argc, char *argv[])
{
  char *ret;
  int num_args_processed = 0;

  if (argc < 3)
  {
    printf ("Invalid arguments for mode 'surface'!\n");
    print_help ();
    exit (EXIT_FAILURE);
  }

  CONFIG.arg_tau_surface = strtod (argv[2], &ret);
  if (*ret != '\0')
  {
    printf ("Unable to convert '%s' into a number for mode 'surface'\n", argv[2]);
    exit (EXIT_FAILURE);
  }
  if (CONFIG.arg_tau_surface <= 0)
  {
    printf ("Invalid value for optical depth, must be positive and non-zero\n");
    exit (EXIT_FAILURE);
  }

  num_args_processed += 1;

  return num_args_processed;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]  argc  The number of arguments provided
 * @param[in]  argv  The command line arguments
 *
 * @details
 *
 * ************************************************************************** */

static int
parse_shared_args (int argc, char *argv[])
{
  char *ret;
  int num_args_processed = 0;

  // Start from 2 to ignore program name and run mode choice
  for (int i = 2; i < argc; ++i)
  {
    if (strcmp (argv[i], "-") == 0)
    {
      printf ("Invalid argument %s\n", argv[i]);
      print_help ();
      exit (EXIT_FAILURE);
    }
    else if (!strcmp (argv[i], "--smax"))
    {
      SMAX_FRAC = strtod (argv[i + 1], &ret);
      if (*ret != '\0')
      {
        printf ("Unable to convert argument for -smax into a double");
        exit (EXIT_FAILURE);
      }
      i += 1;                   // avoid trying to process the argument for this argument
      num_args_processed += 2;
    }
    else if (!strcmp (argv[i], "--freq-min"))   // NOTE: lower frequency boundary for optical depth spectrum
    {
      CONFIG.arg_freq_min = strtod (argv[i + 1], &ret);
      if (*ret != '\0')
      {
        printf ("Unable to convert argument provided for -freq_min to a double\n");
        exit (EXIT_FAILURE);
      }
      i += 1;
      num_args_processed += 2;
    }
    else if (!strcmp (argv[i], "--freq-max"))   // NOTE: upper frequency boundary for optical depth spectrum
    {
      CONFIG.arg_freq_max = strtod (argv[i + 1], &ret);
      if (*ret != '\0')
      {
        printf ("Unable to convert argument provided for -freq_max to a double\n");
        exit (EXIT_FAILURE);
      }
      i += 1;
      num_args_processed += 2;
    }
    else if (!strcmp (argv[i], "--no-es"))
    {
      CONFIG.ignore_electron_scattering = true;
      num_args_processed += 1;
    }
  }

  return num_args_processed;
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
parse_optd_arguments (int argc, char *argv[])
{
  if (argc < 3)
  {
    printf ("Invalid choice of arguments!\n\n");
    print_help ();
    exit (EXIT_FAILURE);
  }

  // deal with one and done arguments
  if (strcmp (argv[1], "--help") == 0 || strcmp (argv[1], "-h") == 0)
  {
    print_help ();
    exit (EXIT_SUCCESS);
  }
  if (!strcmp (argv[1], "--version"))
  {
    printf ("Python version %s\n", VERSION);
    printf ("Built from git commit %s\n", GIT_COMMIT_HASH);
    if (GIT_DIFF_STATUS)
    {
      printf ("This version was compiled with %d files with uncommited changes\n", GIT_DIFF_STATUS);
    }
    exit (EXIT_SUCCESS);
  }

  // Parse the run mode of the program, and then deal with its arguments
  int num_args_processed = 0;

  if (strcmp (argv[1], "pi") == 0)
  {
    CONFIG.mode = MODE_PHOTOION_EDGES;
    num_args_processed += parse_args_for_mode_pi (argc, argv);
  }
  else if (strcmp (argv[1], "wind") == 0)
  {
    CONFIG.mode = MODE_WIND_SPECTRUM;
    num_args_processed += parse_args_for_mode_wind (argc, argv);
  }
  else if (strcmp (argv[1], "cell") == 0)
  {
    CONFIG.mode = MODE_CELL_SPECTRUM;
    num_args_processed += parse_args_for_mode_cell (argc, argv);
  }
  else if (strcmp (argv[1], "surface") == 0)
  {
    CONFIG.mode = MODE_FIND_SURFACE;
    num_args_processed += parse_args_for_mode_surface (argc, argv);
  }
  else
  {
    printf ("Unknown run mode: %s\n", argv[1]);
    print_help ();
    exit (EXIT_FAILURE);
  }

  // Now we've gotten the run mode set up, we can parse arguments which are
  // shared between modes
  num_args_processed += parse_shared_args (argc, argv);

  // assume the final argument is the root name of the simulation being
  // post-processed
  if (num_args_processed + 2 == argc)
  {
    printf ("All command line arguments have been processed and did not find a root name\n");
    exit (EXIT_FAILURE);
  }
  get_root (files.root, argv[argc - 1]);
}

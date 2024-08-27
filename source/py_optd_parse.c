/** ************************************************************************* */
/**
 * @file  py_optd_parse.c
 * @author Edward Parkinson
 * @date August 2024
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
 * @brief  Help text for the program.
 *
 * ************************************************************************** */

static void
print_help (void)
{
  char *help = "A utility program to analyse the optical depth in a Python model.\n\n"
    "usage: py_optical_depth mode [mode-args] root\n\n" "Here is some more text, TODO this need to be finished\n";
  printf ("%s", help);
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
parse_spectrum_args (int argc, char *argv[], struct CommandlineArguments *arguments)
{
  int num_args_processed = 0;

  for (int i = 2; i < argc; ++i)
  {
    if (!strcmp (argv[i], "--inc"))
    {
      char *stdtod_check;
      int num_angles = 0;

      for (int j = 1; j < MAX_CUSTOM_ANGLES + 2; j++)
      {
        double inclination = strtod (argv[i + j], &stdtod_check);
        if (inclination < 0 || inclination > 90)
        {
          printf ("Invalid inclination angle of %f. Has to be 0 < inclination < 90\n", inclination);
          exit (EXIT_FAILURE);
        }
        if (j > MAX_CUSTOM_ANGLES)
        {
          printf ("Maximum number of inclination angles provided: only %d are allowed\n", MAX_CUSTOM_ANGLES);
          exit (EXIT_FAILURE);
        }
        if (strcmp (argv[i + j], stdtod_check) == 0)    // conversion error indicates no more numbers
        {
          break;
        }
        arguments->inclinations[j - 1] = inclination;
        num_angles += 1;
      }

      if (num_angles == 0)
      {
        printf ("Invalid input: no inclination were provided for --inc option\n");
        exit (EXIT_FAILURE);
      }

      num_args_processed = num_angles + 1;
      i += num_args_processed;
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
parse_cell_spectrum_args (int argc, char *argv[], struct CommandlineArguments *arguments)
{
  int num_args_processed = 0;

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
parse_surface_args (int argc, char *argv[], struct CommandlineArguments *arguments)
{
  char *check = NULL;
  int num_args_processed = 1;

  (void) argc;
  (void) arguments;

  if (argc < 3)
  {
    printf ("Invalid number of arguments for mode 'surface'\n");
    print_help ();
    exit (EXIT_FAILURE);
  }

  TAU_DEPTH = strtod (argv[2], &check);
  if (*check != '\0')
  {
    printf ("Unable to convert '%s' into a number for mode 'surface'\n", argv[2]);
    exit (EXIT_FAILURE);
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
parse_photoion_args (int argc, char *argv[], struct CommandlineArguments *arguments)
{
  int num_args_processed = 0;

  (void) argc;
  (void) arguments;

  if (!strcmp (argv[2], "--ion"))       //NOTE: extract column density for specific ion, rather than H
  {
    char *check;
    COLUMN_MODE = COLUMN_MODE_ION;
    COLUMN_MODE_ION_NUMBER = (int) strtol (argv[3], &check, 10);
    if (*check != '\0')
    {
      printf ("Unable to convert argument provided for --ion into an integer\n");
      exit (EXIT_FAILURE);
    }
    if (COLUMN_MODE_ION_NUMBER < 0)
    {
      printf ("Argument for --ion cannot be negative\n");
      exit (EXIT_FAILURE);
    }
    num_args_processed = 2;
  }

  if (COLUMN_MODE == COLUMN_MODE_ION)
  {
    if (COLUMN_MODE_ION_NUMBER < 0 || COLUMN_MODE_ION_NUMBER > nions - 1)
    {
      printf ("The ion number %i is an invalid ion number: there are %i ions which have been loaded\n", COLUMN_MODE_ION_NUMBER, nions);
      exit (EXIT_FAILURE);
    }

    printf ("Extracting column density for %s %i\n", ele[ion[COLUMN_MODE_ION_NUMBER].nelem].name, ion[COLUMN_MODE_ION_NUMBER].istate);
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
parse_global_args (int argc, char *argv[], struct CommandlineArguments *arguments)
{
  char *check = NULL;
  int num_args_processed = 0;

  for (int i = 2; i < argc; i++)
  {
    if (!strncmp (argv[i], "-", 1))
    {
      printf ("Unknown argument %s\n", argv[i]);
      print_help ();
      exit (EXIT_FAILURE);
    }
    else if (!strcmp (argv[i], "--smax"))
    {
      SMAX_FRAC = strtod (argv[i + 1], &check);
      if (*check != '\0')
      {
        printf ("Unable to convert argument for -smax into a double");
        exit (EXIT_FAILURE);
      }
      num_args_processed = i++;
    }
    else if (!strcmp (argv[i], "--domain"))     //NOTE: change the lanching domain for photons
    {
      N_DOMAIN = (int) strtol (argv[i + 1], &check, 10);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for --domain into an integer\n");
        exit (EXIT_FAILURE);
      }
      num_args_processed = i++;
    }
    else if (!strcmp (argv[i], "--freq-min"))   //NOTE: lower frequency boundary for optical depth spectrum
    {
      arguments->freq_min = strtod (argv[i + 1], &check);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -freq_min to a double\n");
        exit (EXIT_FAILURE);
      }
      num_args_processed = i++;
    }
    else if (!strcmp (argv[i], "--freq-max"))   //NOTE: upper frequency boundary for optical depth spectrum
    {
      arguments->freq_max = strtod (argv[i + 1], &check);
      if (*check != '\0')
      {
        printf ("Unable to convert argument provided for -freq_max to a double\n");
        exit (EXIT_FAILURE);
      }
      num_args_processed = i++;
    }
    else if (!strcmp (argv[i], "--no-es"))
    {
      RUN_MODE = MODE_IGNORE_ELECTRON_SCATTERING;
      num_args_processed = i++;
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

extern struct CommandlineArguments
get_arguments (int argc, char *argv[])
{
  struct CommandlineArguments arguments;
  arguments.freq_min = arguments.freq_max = 0;
  arguments.n_freq = NUM_FREQUENCY_BINS;
  for (int i = 0; i < MAX_CUSTOM_ANGLES; ++i)
  {
    arguments.inclinations[i] = -1.0;
  }

  if (argc < 3)
  {
    printf ("Invalid number of arguments\n");
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
      printf ("This version was compiled with %d files with uncommited changes\n", GIT_DIFF_STATUS);
    exit (EXIT_SUCCESS);
  }

  // deal with run mode and their arguments
  int num_args_processed = 0;
  char mode_string[LINELENGTH];
  strncpy (mode_string, argv[1], LINELENGTH);

  if (strcmp (mode_string, "edges") == 0)
  {
    RUN_MODE = MODE_PHOTOION;
    num_args_processed += parse_photoion_args (argc, argv, &arguments);
  }
  else if (strcmp (mode_string, "wind") == 0)
  {
    RUN_MODE = MODE_SPECTRUM;
    num_args_processed += parse_spectrum_args (argc, argv, &arguments);
  }
  else if (strcmp (mode_string, "cell") == 0)
  {
    RUN_MODE = MODE_CELL_SPECTRUM;
    num_args += parse_cell_spectrum_args (argc, argv, &arguments);
  }
  else if (strcmp (mode_string, "surface") == 0)
  {
    RUN_MODE = MODE_SURFACE;
    num_args_processed += parse_surface_args (argc, argv, &arguments);
  }
  else
  {
    printf ("Unknown mode: %s\n", mode_string);
    print_help ();
    exit (EXIT_FAILURE);
  }

  // deal with global arguments which have multiple uses
  num_args_processed += parse_global_args (argc, argv, &arguments);

  if (num_args_processed + 2 == argc)
  {
    printf ("All command line arguments have been processed and did not find a root name\n");
    exit (EXIT_FAILURE);
  }

  // assume the final argument is the root name of the simulation
  get_root (files.root, argv[argc - 1]);

  return arguments;
}

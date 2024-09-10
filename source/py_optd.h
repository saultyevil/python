/** ************************************************************************* */
/**
 * @file     py_optd.h
 * @author   Edward Parkinson
 * @date     February 2021
 *
 * @details
 *
 * This header file contains the main constants, macros and types used in
 * py_optd.
 *
 * ************************************************************************** */

#ifndef PY_OPTD_H
#define PY_OPTD_H

#include <stdbool.h>

#define MAX_CUSTOM_ANGLES 64
#define MAX_ANGLES 128
#define NUM_FREQUENCY_BINS 10000
#define NAMELEN 32

// Error message macro, adds the file name and line to the start of the error
// message

#define print_error(fmt, ...)                         \
  {                                                   \
    fprintf(stderr, "(%s:%i): ", __FILE__, __LINE__); \
    fprintf(stderr, fmt, ##__VA_ARGS__);              \
  }

/** Structure to hold the angles to extract the optical depth/column density from
 */

typedef struct SightLines
{
  char name[NAMELEN];
  double angle;
  double direction_vector[3];
} SightLines_t;

/** Structure to hold the name and frequency of a photoionization edge to
 *evaluate the optical depth at
 */

typedef struct PIEdge
{
  char name[50];
  double freq;
} PIEdge_t;

/** Structure to save the angle and position of the electron scattering
 * photosphere surface
 */

typedef struct Pos
{
  double angle;
  double x, y, z;
} Pos_t;

/** Enumerator used to control the column density which is extracted, i.e. by
 * default mass density/N_H is extracted by the density of an ion can also
 * be extracted
 */

enum ColumnDensityEnum
{
  COLUMN_MODE_RHO,
  COLUMN_MODE_ION,
};

enum RunModeEnum
{
  MODE_PHOTOION_EDGES,
  MODE_WIND_SPECTRUM,
  MODE_CELL_SPECTRUM,
  MODE_FIND_SURFACE,
};

struct Config
{
  enum RunModeEnum mode;
  enum ColumnDensityEnum column_density;
  bool ignore_electron_scattering;
  int column_density_ion_number;
  int domain;
  double arg_freq_min;
  double arg_freq_max;
  int arg_num_inc;
  double arg_inclinations[MAX_ANGLES];
  double arg_tau_surface;
  int arg_wind_elem;
  int arg_wind_i;
  int arg_wind_j;
} CONFIG;

// Function prototypes

void parse_optd_arguments(int argc, char *argv[]);
int initialize_inclination_angles(struct SightLines *inclinations);
int initialise_photon_packet(PhotPtr photon, double frequency, double *direction);
int integrate_tau_across_wind(PhotPtr photon, double *c_column_density, double *c_optical_depth);
int integrate_tau_across_cell(PhotPtr photon, double *column_density_out, double *optical_depth_out);
int print_optical_depths(SightLines_t *inclinations, int n_inclinations, PIEdge_t edges[], int n_edges, double *optical_depth,
                         double *column_density);
int save_optical_depth_spectrum(SightLines_t *inclinations, int n_inclinations, double *tau_spectrum, double freq_min, double d_freq);
void write_photosphere_location_to_file(Pos_t *positions, int n_angles);
int initialize_wind_structures(void);
void set_frequency_range(void);
int calculate_photoionization_optical_depths(void);
int create_optical_depth_spectrum(void);
int find_optical_depth_surface(void);

#endif

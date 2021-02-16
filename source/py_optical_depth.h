
/** ************************************************************************* */
/**
 * @file     py_optical_depth.h
 * @author   Edward Parkinson
 * @date     February 2021
 *
 * @brief    Globals used during optical depth diagnostics.
 *
 * @details
 *
 * Do not include this in any file other than py_optical_depth_sub.c, otherwise
 * there will be double definition errors. If you want to include it in
 * multiple, place MAXDIFF and N_PI_EDGES into py_optical_depth_sub.c to avoid
 * the errors. The reason they are in this header file, and why it exists, is
 * to keep things a bit cleaner in py_optical_depth_sub.c
 *
 * ************************************************************************** */

#define DOMAIN_TO_CONSIDER 0  // For now, we only care about photons starting in domain 0
#define N_FREQ_BINS 25000

// Inclination angles structure

typedef struct SightLines_s
{
  char name[50];
  double lmn[3];
} SightLines_t;

SightLines_t *INCLINATION_ANGLES;
int N_INCLINATION_ANGLES;         // The number of inclination angles

// PI edges structure

typedef struct PIEdges_s
{
  char name[50];
  double freq;
} PIEdges_t;

#define N_PI_EDGES (int) (sizeof PI_EDGES_TO_MEASURE / sizeof *PI_EDGES_TO_MEASURE)  // The number of optical depths for the simple calculation

// Control the column density extracted

enum COLUMN_DENSITY
{
  COLUMN_MODE_RHO,
  COLUMN_MODE_ION,
};

int COLUMN_MODE;
int COLUMN_MODE_ION_NUMBER;

/** ************************************************************************* */
/**
 * @file     tau_diag.h
 * @author   Edward Parkinson
 * @date     May 2019
 *
 * @brief    The main resting place for all constants and global vairables
 *           related to providing optical depth diagnostic.
 *
 * @details
 *
 *
 * ************************************************************************** */

#ifndef PYTHON_TAU_DIAG_H
#define PYTHON_TAU_DIAG_H

/*
 * SPEC_OPAC
 * ---------
 * Enumerator to label any special opacities which are used. This is mainly for
 * the use of the Rosseland and Planck mean opacities.
 */

enum SPEC_OPAC
{
  NORMAL_TAU = 1,
  ROSSEL_MEAN = 2,
  PLANCK_MEAN = 3,
};

/*
 * TAU_DIAG_OPACS
 * --------------
 * This is the global variable type used to track the name of opacities for the
 * optical depth diagnostics. This includes the frequency of each photon, nu,
 * and the name of the opacity.
 */

typedef struct tau_diag_opacities
{
  double nu;
  char name[50];
} Opacities;

Opacities TAU_DIAG_OPACS[] = {
  {-PLANCK_MEAN, "Planck"},
  {-ROSSEL_MEAN, "Rosseland"},
  {3.387485e+15, "Lyma_885A"},
  {8.293014e+14, "Bal_3615A"},
  {1.394384e+16, "HeII_215A"}
};

int const N_TAU = sizeof TAU_DIAG_OPACS / sizeof *TAU_DIAG_OPACS;

/*
 * tau_diag_observers
 * ------------------
 * This is the global variable type used to track the name and cosine direction
 * vectors of the inclinations of which observers are placed. Hence, the members
 * of this type will be the various inclination angles for integrated optical
 * depths.
 */

typedef struct observer_angles
{
  char name[50];
  double lmn[3];
} Observers;

Observers *tau_diag_observers;

int N_ANGLES;  // The number of inclination angles
double print_xloc;  // Used as a lazy way to pass the vertical photon's x location
#define DEFAULT_DOMAIN 0  // For now, assume we only care about domain 0 diagnostics

#endif

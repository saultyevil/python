/** ************************************************************************* */
/**
 * @file     tau_diag.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    The main resting place for most functions related to providing
 *           optical depth diagnostics.
 *
 * @details
 *
 * ### HOW TO ADD MORE OPTICAL DEPTHS ###
 *
 * In theory it should be relatively straight forwards to add more optical
 * depths if desired. One should only need to update the the TAU_DIAG_OPACS
 * array in tau_diag.h.
 *
 * All one should be required to do is to add a name for the optical depth and
 * any relevant frequency for photons to take. * This is done as the function
 * radiation is called to calculate the opacity for a photon of a given
 * frequency in the current cell.
 *
 * If one wishes to add a "special" opacity, which is an opacity which will
 * require it's own function to calculate, such as the Rosseland mean opacity,
 * then one should update the enumerator SPEC_OPAC in tau_diag.h with a new
 * label for that opacity, as well as the frequency to be -ENUMERATOR in the
 * TAU_DIAG_OPACS array. The frequency is set such as way so as in the function
 * tau_diag, one can update the if statement to include the new opac type.
 *
 * If none of this makes sense, just try to follow how everything else is done!
 *
 * ### Programming Notes ###
 *
 * The macro __func__ is used in here for error reporting usage, but it is only
 * defined since the C99 standard... but what if people use C90? oh no! Not
 * sure if I should really care?
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"

/*
 * tau_diag_observers
 * ------------------
 * This is the global variable type used to track the name and cosine direction
 * vectors of the inclinations of which observers are placed. Hence, the members
 * of this type will be the various inclination angles for integrated optical
 * depths.
 */

struct observer_angles
{
  char name[50];
  double lmn[3];
} *tau_diag_observers;

/*
 * SPEC_OPAC
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
  {3.387485e+15, "LymanEdge_885A"},
  {8.293014e+14, "BalmerEdge_3615A"},
  {1.394384e+16, "HeIIEdge_215A"}
};

int const N_TAU = sizeof TAU_DIAG_OPACS / sizeof *TAU_DIAG_OPACS;

int N_ANGLES;                   // The number of inclination angles
#define DEFAULT_DOMAIN 0        // For now, assume we only care about domain 0 diagnostics

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
 * Historically, this used to also print to file. But I found that it was mostly
 * useless to do this, considering the information is already in the diag files
 * by using Log ().
 *
 * ************************************************************************** */

void
print_tau_angles (double *tau_store, double *col_den_store)
{
  int itau;
  int ispec;
  int line_len;
  int const MAX_COL = 120;
  double tau;
  double col_den;
  char tmp_str[LINELENGTH];
  char observer_name[LINELENGTH];

  if (rank_global != 0)
    return;

  Log ("Optical depths along the defined line of sights:\n\n");

  for (ispec = 0; ispec < N_ANGLES; ispec++)
  {
    strcpy (observer_name, tau_diag_observers[ispec].name);

    Log ("%s: ", observer_name);

    col_den = col_den_store[ispec];
    Log ("Column mass density: %3.2e g/cm^-2\n", col_den);

    line_len = 0;
    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = tau_store[ispec * N_TAU + itau];
      line_len += sprintf (tmp_str, "tau_%-9s: %3.2e  ", TAU_DIAG_OPACS[itau].name, tau);

      if (line_len > MAX_COL)
      {
        line_len = 0;
        Log ("\n");
      }

      Log ("%s", tmp_str);
    }

    Log ("\n\n");
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
write_tau_spectrum_to_file (double *tau_spectrum, double freq_min, double dfreq)
{
  int ispec;
  int ifreq;
  double wavelength;
  double frequency;
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
   * Write the file header - nothing fancy
   */

  fprintf (tau_spec_file, "# Lambda ");
  for (ispec = 0; ispec < N_ANGLES; ispec++)
    fprintf (tau_spec_file, "%s ", tau_diag_observers[ispec].name);
  fprintf (tau_spec_file, "\n");

  /*
   * Write out the tau spectrum for each inclination angle
   */

  frequency = freq_min;

  for (ifreq = 0; ifreq < NWAVE; ifreq++)
  {
    wavelength = (VLIGHT / frequency) / ANGSTROM;
    fprintf (tau_spec_file, "%e ", wavelength);

    for (ispec = 0; ispec < N_ANGLES; ispec++)
      fprintf (tau_spec_file, "%e ", tau_spectrum[ispec * NWAVE + ifreq]);

    fprintf (tau_spec_file, "\n");

    frequency += dfreq;
  }

  if (fclose (tau_spec_file))
  {
    Error ("%s:%s:%i: could not close tau spectrum diag file\n", __FILE__, __func__, __LINE__);
    Exit (1);
  }
}

/* ************************************************************************* */
/**
 * @brief           Calculate the total optical depth a photon experiences
 *                  across the cell of distance smax.
 *
 * @param[in]       WindPtr   w           A pointer to the entire wind
 * @param[in]       PhotPtr   pextract    The photon packet to extract
 * @param[in]       int       opac_type   The opacity to use to calculate the
 *                                        optical depth
 * @param[in,out]   double    *col_den    The column density the photon has
 *                                        translated through
 * @param[in,out]   double    *tau        The optical depth experienced by the
 *                                        photon
 *
 * @return          int       istat       The current photon status or
 *                                        EXIT_FAILURE on failure.
 *
 * @details
 *
 * This function is concerned with finding the opacity of the photon's current
 * cell, the distance the photon can move in the cell and hence it increments
 * the optical depth tau a photon has experienced as it moves through the wind.
 *
 * The opacity is usually calculated using the function radiation which the
 * extract photon is passed to, as well as the distance it can move, smax. In
 * special cases, such as when a mean opacity is being used, radiation is not
 * used to find the opacity but specialised functions are called instead.
 *
 * The photon status istat is updated and returned for some reason.
 *
 * TODO: implement mean opacities
 * TODO: replace mass column density with electron or H densities
 *
 * ************************************************************************** */

int
calculate_tau (WindPtr w, PhotPtr pextract, int opac_type, double *col_den, double *tau)
{
  int istat;
  int nplasma;
  double density;
  double smax;
  double kappa_tot;

  PlasmaPtr plasma_cell;

  if ((pextract->grid = where_in_grid (w[pextract->grid].ndom, pextract->x)) < 0)
  {
    Error ("%s:%s:%i: pextract is not in grid\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  /*
   * Determine where in the plasma grid the cell is. This is required so a
   * column density can be calculated
   */

  nplasma = w[pextract->grid].nplasma;
  plasma_cell = &plasmamain[nplasma];
  density = plasma_cell->rho;

  /*
   * smax should be the transverse distance of the cell the photon is in
   */

  if ((smax = find_smax (pextract)) < 0)
  {
    Error ("%s:%s:%i: abnormal value of smax for photon\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  /*
   * Depending on the type of opacity in use, kappa_tot is updated with the
   * opacity of the cell which is the total opacity of that cell. Note that
   * radiation is the opacity of the cell in non-macro atom mode.
   */

  if (opac_type == ROSSEL_MEAN)
  {
    return P_ESCAPE;
  }
  else if (opac_type == PLANCK_MEAN)
  {
    return P_ESCAPE;
  }
  else
  {
    kappa_tot = radiation (pextract, smax);
  }

  /*
   * Increment the optical depth and column density variables and move the
   * photon to the edge of the cell
   */

  *col_den += smax * density;
  *tau += smax * kappa_tot;
  move_phot (pextract, smax);

  istat = pextract->istat;

  return istat;
}

/* ************************************************************************* */
/**
 * @brief           Extract the optical depth the photon packet porig must
 *                  travel through to reach the observer.
 *
 * @param[in]       WindPtr   w          A pointer to the entire wind
 * @param[in]       PhotPtr   porig      The photon packet to extract
 * @param[in]       int       opac_type  An indicator of wheter to use a mean
 *                                       opacity or not
 * @param[out]      double    *col_den   The column depth of the extracted
 *                                       photon angle
 * @param[out]      double    *tau       The optical depth from photon origin to
 *                                       the observer
 *
 * @return          int                  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * The direction of the observer is set in porig.lmn and should be set prior
 * to passing a photon to this function. If any values of lmn are negative, the
 * photon will translate in the negative direction. However, have no fear as this
 * is normal and is fine due to the assumed symmetry of models in Python.
 *
 * ************************************************************************** */

int
extract_tau (WindPtr w, PhotPtr porig, int opac_type, double *col_den, double *tau)
{
  int rc;
  int ndom;
  int istat;
  int n_trans_space;
  int const max_trans_space = 10;

  double norm[3];

  struct photon pextract;

  // TODO: remove when mean opacities are implemeneted
  if (opac_type == ROSSEL_MEAN || opac_type == PLANCK_MEAN)
    return EXIT_FAILURE;

  istat = P_INWIND;             // assume photon is in wind for initialisation reasons
  stuff_phot (porig, &pextract);

  n_trans_space = 0;

  while (istat == P_INWIND)
  {
    if (where_in_wind (pextract.x, &ndom) < 0)  // Move the photon in the grid
    {
      translate_in_space (&pextract);
      if (++n_trans_space > max_trans_space)
      {
        Error ("%s:%s:%i: photon transport ended due to too many translate_in_space\n", __FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
      }
    }
    else if ((pextract.grid = where_in_grid (ndom, pextract.x)) >= 0)   // Move the photon in the wind
    {
      rc = calculate_tau (w, &pextract, opac_type, col_den, tau);
      if (rc)
      {
        pextract.istat = -1;
        return EXIT_FAILURE;
      }
    }
    else
    {
      Error ("%s:%s:%i: photon in unknown location, grid stat %i\n", __FILE__, __func__, __LINE__, pextract.grid);
      pextract.istat = -1;
      return EXIT_FAILURE;
    }

    /*
     * Now that the photon has been translated, update istat using the walls.
     * This should return istat != P_INWIND when the photon has escaped or smt
     */

    istat = walls (&pextract, porig, norm);
  }

  if (istat == P_HIT_STAR || istat == P_HIT_DISK)
  {
    Error ("%s:%s:%i: photon hit star or disk, istat = %i\n", __FILE__, __func__, __LINE__, istat);
    return EXIT_FAILURE;
  }

  if (walls (&pextract, porig, norm) == -1)
  {
    Error ("%s:%s:%i: abnormal return from walls for photon\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief           Reposition a photon depending on various criteria
 *
 * @param[in, out]  PhotPtr pout    The photon to reposition
 *
 * @return          void
 *
 * @details
 *
 * This function assumes that we only care about DEFAULT_DOMAIN and whatever
 * that has been set as.
 *
 * ************************************************************************** */

void
reposition_tau_photon (PhotPtr pout)
{
  int icell;
  int wind_index;
  int plasma_index;
  int static xloc_init = 0;

  double x_loc;
  double static upward_x_loc;
  double plasma_density;
  double max_plasma_density;

  /*
   * Hacky code for when photon is pointing upwards - this photon will not
   * be launched from the origin, but will be launched from the most dense
   * cell at the bottom of the disk.
   */

  if (xloc_init == 0)
  {
    max_plasma_density = 0;
    for (icell = 0; icell < zdom[DEFAULT_DOMAIN].ndim - 1; icell++)
    {
      wind_index = icell * zdom[DEFAULT_DOMAIN].mdim + zdom[DEFAULT_DOMAIN].mdim;
      x_loc = wmain[wind_index].x[0];
      plasma_index = wmain[wind_index].nplasma;
      plasma_density = plasmamain[plasma_index].rho;

      if (plasma_density > max_plasma_density)
      {
        upward_x_loc = x_loc + wmain[wind_index].dfudge;
        max_plasma_density = plasma_density;
      }
    }

    xloc_init = 1;
  }

  if (pout->lmn[0] == 0 && pout->lmn[1] == 0 && pout->lmn[2] == 1)
    pout->x[0] = upward_x_loc;

  /*
   * When the photon is pointing in the negative x-direction, i.e. when
   * ptau.lmn[0] < 0, then it's very likely that the photon will hit the
   * central object and also not travel through the base of this wind. In
   * this case, the photon is simply placed on the negative side of the
   * x-axis which is absolutely fine due to the rotational symmetry in
   * Python
   */

  if (pout->lmn[0] < 0)
    pout->x[0] *= -1;
}

/* ************************************************************************* */
/**
 * @brief           Generate a photon packet with a given frequency nu for use
 *                  with the optical depth diagnostic routines.
 *
 * @param[in,out]   PhotPtr  pout     The photon packet to initialise
 *
 * @param[in]       double   nu       The frequency of the photon packet
 *
 * @return          int               EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * This routine assumes that the user wishes for the photon to be generated from
 * the inner-most radius of the disk or from the central source of the
 * simulation. Note that photons are initialised with a weight of f_tot as
 * photons are required to have weight, but since the diagnostic functions do
 * not care about the weight of the photon, it is set to something large.
 *
 * ************************************************************************** */

int
create_tau_diag_phot (PhotPtr pout, double nu)
{
  if (nu < 0)
  {
    Error ("%s:%s:%i: photon can't be created with negative nu\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  pout->freq = pout->freq_orig = nu;
  pout->origin = pout->origin_orig = PTYPE_DISK;
  pout->istat = P_INWIND;
  pout->w = pout->w_orig = geo.f_tot;
  pout->x[0] = 1.1 * geo.rstar;
  pout->x[1] = 0.0;
  pout->x[2] = 0.0;
  pout->tau = 0.0;

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief           Initialise the viewing angles for the optical depth
 *                  diagnostic routines.
 *
 * @param           void
 *
 * @return          void
 *
 * @details
 *
 * This purpose of this function is to initialise the angles of which to extract
 * optical depth information from. If xxspec is initialised, then the optical
 * depth information will be extracted from the observer angles provided in the
 * parameter file for spectrum generation. Otherwise, the routines will instead
 * use a set of pre-defined angles.
 *
 * ************************************************************************** */

void
init_tau_diag_angles (void)
{
  int iangle;
  int memory_req;

  double default_phase = 0.5;
  double default_angles[] = { 0.0, 10.0, 40.0, 60.0, 80.0, 90.0 };
  int const n_default_angles = sizeof default_angles / sizeof default_angles[0];

  /*
   * Use the angles specified for by the user for spectrum generation, this
   * requires for xxspec to be initialised.
   */

  if (geo.nangles && xxspec != NULL)
  {
    Log_silent ("Using inclinations provided for spectrum cycles as defined in xxspec\n");

    N_ANGLES = geo.nangles;
    tau_diag_observers = calloc (geo.nangles, sizeof *tau_diag_observers);

    if (tau_diag_observers == NULL)
    {
      memory_req = geo.nangles * sizeof *tau_diag_observers;
      Error ("%s:%s:%i: cannot allocate %d bytes for observers array\n", __FILE__, __func__, __LINE__, memory_req);
      Exit (1);
    }
    else
    {
      for (iangle = MSPEC; iangle < MSPEC + geo.nangles; iangle++)
      {
        strcpy (tau_diag_observers[iangle - MSPEC].name, xxspec[iangle].name);
        stuff_v (xxspec[iangle].lmn, tau_diag_observers[iangle - MSPEC].lmn);
      }
    }
  }

  /*
   * If no spectrum cycles are being run, or if the model is being restarted
   * without initialising xxspec, then a default set of angles is used instead.
   */

  else
  {
    Log_silent ("As there are no spectrum cycles or observers defined, a set of default angles will be used instead\n");

    N_ANGLES = n_default_angles;
    tau_diag_observers = calloc (n_default_angles, sizeof *tau_diag_observers);

    if (tau_diag_observers == NULL)
    {
      memory_req = n_default_angles * sizeof *tau_diag_observers;
      Error ("%s:%s:%i: cannot allocate %d bytes for observers array\n", __FILE__, __func__, __LINE__, memory_req);
      Exit (1);
    }
    else
    {
      for (iangle = 0; iangle < n_default_angles; iangle++)
      {
        sprintf (tau_diag_observers[iangle].name, "A%02.0fP%04.2f", default_angles[iangle], default_phase);
        tau_diag_observers[iangle].lmn[0] = sin (default_angles[iangle] / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
        tau_diag_observers[iangle].lmn[1] = sin (default_angles[iangle] / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
        tau_diag_observers[iangle].lmn[2] = cos (default_angles[iangle] / RADIAN);
      }
    }
  }
}

/* ************************************************************************* */
/**
 * @brief           Create spectra of tau vs lambda for each observer angle
 *
 * @param[in]       WindPtr    w         A pointer to the wind
 *
 * @return          void
 *
 * @details
 *
 * This is the main function which will generate the optical depth spectra for
 * each observer angle in xxspec. The algorithm is similar to extract and the
 * tau_diag algorithm which this function is called in.
 *
 * A photon is generated at the central source of the model and is extracted
 * from this location towards the observer where it escapes, where extract_tau
 * returns the integrated optical depth along its path to escape. This is done
 * for a range of photon frequencies to find how optical depth changes with
 * frequency.
 *
 * This processes can take some time compared to tau_diag. But since only NWAVE
 * photons are being generated for each spectrum and the fact that these photons
 * do not interact, this spectrum does not actually take that long.
 *
 * ************************************************************************** */

void
create_tau_spectrum (WindPtr w)
{
  int ispec;
  int ifreq;
  int memory_req;

  double tau;
  double col_den;
  double freq;
  double dfreq;
  double freq_min;
  double freq_max;
  double wave_min;
  double wave_max;
  double *current_observer;
  double *tau_spectrum;

  struct photon ptau;

  Log ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (N_ANGLES * NWAVE, sizeof *tau_spectrum);
  if (!tau_spectrum)
  {
    memory_req = N_ANGLES * NWAVE * sizeof *tau_spectrum;
    Error ("%s:%s:%i: cannot allocate %d bytes for tau_spectrum\n", __FILE__, __func__, __LINE__, memory_req);
    Exit (1);
  }

  /*
   * Define the limits of the spectra in both wavelength and frequency. The
   * size of the frequency and wavelength bins is also defined. Note that the
   * wavelength limits MUST be in units of Angstroms.
   */

  wave_min = 100;
  wave_max = 10000;

  freq_max = VLIGHT / (wave_min * ANGSTROM);
  freq_min = VLIGHT / (wave_max * ANGSTROM);
  dfreq = (freq_max - freq_min) / NWAVE;

  Log_silent ("Minimum wavelength : %.0f Angstroms\n", wave_min);
  Log_silent ("Maximum wavelength : %.0f Angstroms\n", wave_max);
  Log_silent ("Wavelength bins    : %i\n\n", NWAVE);

  /*
   * Now create the optical depth spectra for each observer angle
   */

  for (ispec = 0; ispec < N_ANGLES; ispec++)
  {
    Log ("  - Creating spectrum: %s\n", tau_diag_observers[ispec].name);

    freq = freq_min;
    current_observer = tau_diag_observers[ispec].lmn;

    for (ifreq = 0; ifreq < NWAVE; ifreq++)
    {
      tau = 0;
      col_den = 0;

      if (create_tau_diag_phot (&ptau, freq) == EXIT_FAILURE)
      {
        Log ("%s:%s:%i: skipping photon of frequency %e\n", __FILE__, __func__, __LINE__, freq);
        continue;
      }

      /*
       * Point the photon in the direction of the observer and reposition it if
       * required
       */

      stuff_v (current_observer, ptau.lmn);
      reposition_tau_photon (&ptau);

      if (extract_tau (w, &ptau, NORMAL_TAU, &col_den, &tau) == EXIT_FAILURE)
        tau = -1;

      tau_spectrum[ispec * NWAVE + ifreq] = tau;
      freq += dfreq;
    }
  }

  write_tau_spectrum_to_file (tau_spectrum, freq_min, dfreq);
  free (tau_spectrum);
}

/* ************************************************************************* */
/**
 * @brief           The main controlling function for the optical depth
 *                  diagnostics.
 *
 * @param[in]       WindPtr    w       A pointer to the wind
 *
 * @return          void
 *
 * @details
 *
 * This is the main function which will control the procedure for calculating
 * various diagnostic numbers for the optical depth's experienced in the current
 * model. Namely, this function aims to show the total integrated optical depth
 * to each observer angle using the following optical depths:
 *
 *  - Rosseland mean
 *  - Planck mean
 *  - Lymann edge
 *  - Balmer edge
 *  - Helium II edge
 *
 * Once these integrated optical depths have been calculated for each angle, a
 * spectrum of optical depth vs wavelength is created for each angle.
 *
 * The aim of these diagnostic numbers it to provide some sort of metric on the
 * optical thickness of the current model.
 *
 * TODO: method for calling this for individual cells instead of the entire wind
 *
 * ************************************************************************** */

void
tau_integrate_angles (WindPtr w)
{
  int itau;
  int ispec;
  int opac_type;
  int memory_req;

  double nu;
  double tau;
  double col_den;
  double *current_observer;

  double *tau_store;
  double *col_den_store;

  struct photon ptau;

  tau_store = calloc (N_ANGLES * N_TAU, sizeof *tau_store);
  if (!tau_store)
  {
    memory_req = N_ANGLES * N_TAU * sizeof *tau_store;
    Error ("%s:%s:%i: cannot allocate %d bytes for tau_store\n", __FILE__, __func__, __LINE__, memory_req);
    Exit (1);
  }

  col_den_store = calloc (N_ANGLES, sizeof *tau_store);
  if (!col_den_store)
  {
    memory_req = N_ANGLES * sizeof *tau_store;
    Error ("%s:%s:%i: cannot allocate %d bytes for col_den_store\n", __FILE__, __func__, __LINE__, memory_req);
    Exit (1);
  }

  /*
   * Now extract the optical depths and mass column densities
   */

  for (ispec = 0; ispec < N_ANGLES; ispec++)
  {
    current_observer = tau_diag_observers[ispec].lmn;

    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = 0.0;
      col_den = 0.0;
      nu = TAU_DIAG_OPACS[itau].nu;

      /*
       * If the frequency is negative, a mean opacity will be used instead
       * TODO: not yet implemented :-(
       */

      if (nu == -ROSSEL_MEAN)
      {
        opac_type = ROSSEL_MEAN;
        nu = 1e15;
      }
      else if (nu == -PLANCK_MEAN)
      {
        opac_type = PLANCK_MEAN;
        nu = 1e15;
      }
      else
      {
        opac_type = NORMAL_TAU;
      }

      /*
       * Create the tau diag photon and point it towards the extract angle
       */

      if (create_tau_diag_phot (&ptau, nu) == EXIT_FAILURE)
      {
        Log ("%s:%s:%i: skipping photon of frequency %e\n", __FILE__, __func__, __LINE__, nu);
        tau_store[ispec * N_TAU + itau] = -1;
        col_den_store[ispec] = -1;
        continue;
      }

      /*
       * Point the photon in the direction of the observer and reposition it if
       * required
       */

      stuff_v (current_observer, ptau.lmn);
      reposition_tau_photon (&ptau);

      /*
       * Extract tau and add it to the tau storage array
       */

      if (extract_tau (w, &ptau, opac_type, &col_den, &tau) == EXIT_FAILURE)
      {
        col_den = -1;
        tau = -1;
      }

      tau_store[ispec * N_TAU + itau] = tau;
      col_den_store[ispec] = col_den;
    }
  }

  print_tau_angles (tau_store, col_den_store);

  free (tau_store);
  free (col_den_store);
}

/* ************************************************************************* */
/**
 * @brief           Main control function for create optical depth diagnostics.
 *
 * @param           WindPtr w       A pointer to the wind
 *
 * @return          void
 *
 * @details
 *
 * This function is the main steering function for generating the optical depth
 * diagnostics. In MPI mode, all MPI processes should be doing an identical
 * amount of work, however only the master process (0) will write any diagnostics
 * out to file. To avoid any complications due to the processes being out of
 * sync with one another, a MPI_Barrier is used at the end to allow the processes
 * to resynchronise.
 *
 * ************************************************************************** */

void
tau_diag_main (WindPtr w)
{

  xsignal (files.root, "%-20s Optical depth diagnostics beginning\n", "TAU");

  /*
   * Switch to extract mode - this is generally not great practise but it
   * enables radiation to return the opacity in a cell without updating
   * the monte carlo estimators
   */

  geo.ioniz_or_extract = 0;

  init_tau_diag_angles ();
  tau_integrate_angles (w);
  create_tau_spectrum (w);

  geo.ioniz_or_extract = 1;

  if (tau_diag_observers != NULL)
    free (tau_diag_observers);

#ifdef MPI_ON
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  xsignal (files.root, "%-20s Optical depth diagnostics completed\n", "TAU");
  Log (" Completed optical depth diagnostics.  The elapsed TIME was %f\n", timer ());
}

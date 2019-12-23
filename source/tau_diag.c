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
 * OBSERVERS
 * This is the global variable type used to track the name and cosine direction
 * vectors of the inclinations of which observers are placed. Hence, the members
 * of this type will be the various inclination angles for integrated optical
 * depths.
 */

struct observer_angles
{
  char name[50];
  double lmn[3];
} *OBSERVERS;

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

int const N_TAU = sizeof TAU_DIAG_OPACS / sizeof *TAU_DIAG_OPACS;       // The number of optical depths for the simple calculation

int N_ANGLES;                   // The number of inclination angles
#define DEFAULT_DOMAIN 0        // For now, assume we only care about photons starting in domain 0

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
print_tau_angles (const double *tau_store, const double *col_den_store)
{
  int itau, ispec;
  int line_len;
  int const MAX_COL = 120;
  double tau, column;
  char tmp_str[LINELENGTH];
  char observer_name[LINELENGTH];

  if (rank_global != 0)
    return;

  Log ("Optical depths along the defined line of sights:\n\n");

  for (ispec = 0; ispec < N_ANGLES; ispec++)
  {
    strcpy (observer_name, OBSERVERS[ispec].name);

    Log ("%s: ", observer_name);

    column = col_den_store[ispec];
    Log ("Hydrogen column density: %3.2e cm^-2\n", column / MPROT);

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
write_tau_spectrum_to_file (const double *tau_spectrum, double freq_min, double dfreq)
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

  if ((tau_spec_file = fopen (tau_spec_filename, "w")) == NULL)
  {
    Error ("%s:%s:%i: unable to open tau spectrum diag file\n", __FILE__, __func__, __LINE__);
    Exit (1);
  }

  /*
   * Write the file header - nothing fancy
   */

  fprintf (tau_spec_file, "# Lambda Freq. ");
  for (ispec = 0; ispec < N_ANGLES; ispec++)
    fprintf (tau_spec_file, "%s ", OBSERVERS[ispec].name);
  fprintf (tau_spec_file, "\n");

  /*
   * Write out the tau spectrum for each inclination angle
   */

  frequency = freq_min;

  for (ifreq = 0; ifreq < NWAVE; ifreq++)
  {
    wavelength = VLIGHT / frequency / ANGSTROM;
    fprintf (tau_spec_file, "%e %e ", wavelength, frequency);

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
 * ************************************************************************** */

int
calculate_tau (WindPtr w, PhotPtr pextract, double *col_den, double *tau)
{
  int istat;
  int ndom, nplasma;
  double kappa_tot, density;
  double smax;
  double v_inner[3], v_outer[3], v_check[3];
  double v_extract, v_far, vc, vch;
  double freq_inner, freq_outer;
  double mean_freq;

  WindPtr wind_cell;
  PlasmaPtr plasma_cell;
  struct photon p_far, p_mid;

  if ((pextract->grid = where_in_grid (w[pextract->grid].ndom, pextract->x)) < 0)
  {
    Error ("%s:%s:%i: pextract is not in grid\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  /*
   * Determine where in the plasma grid the cell is. This is required so a
   * column density can be calculated
   */

  wind_cell = &w[pextract->grid];
  ndom = wind_cell->ndom;
  nplasma = w[pextract->grid].nplasma;
  plasma_cell = &plasmamain[nplasma];
  density = plasma_cell->rho;

  /*
   * smax should be the transverse distance of the cell the photon is in
   */

  smax = find_smax (pextract);
  if (smax < 0)
  {
    Error ("%s:%s:%i: abnormal value of smax for photon\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  /*
   * Calculate the velocity at both sides of the cell
   */

  stuff_phot (pextract, &p_far);
  move_phot (&p_far, smax);
  vwind_xyz (ndom, pextract, v_inner);
  v_extract = dot (pextract->lmn, v_inner);
  vwind_xyz (ndom, &p_far, v_outer);
  v_far = dot (p_far.lmn, v_outer);

  /*
   * Check to see if the velocity is monotonic across the cell by calculating
   * the velocity at the midpoint of the path. If it is not monotonic, then
   * reduce smax
   */

  vc = VLIGHT;
  while (vc > VCHECK && smax > DFUDGE)
  {
    stuff_phot (pextract, &p_mid);
    move_phot (&p_mid, smax / 2.);
    vwind_xyz (ndom, &p_mid, v_check);
    vch = dot (p_mid.lmn, v_check);
    vc = fabs (vch - 0.5 * (v_extract + v_far));
    if (vc > VCHECK)
    {
      stuff_phot (&p_mid, &p_far);
      smax *= 0.5;
      v_far = vch;
    }
  }

  /*
   * Doppler shift the photons from the global to the local frame of rest
   */

  freq_inner = pextract->freq * (1.0 - v_extract / VLIGHT);
  freq_outer = p_far.freq * (1.0 - v_far / VLIGHT);
  mean_freq = (freq_inner + freq_outer) / 2.0;

  /*
   * Now we can finally calculate the opacity due to all the continuum
   * processes. In macro-atom mode, we need to calculate the continuum opacity
   * using kappa_bf and kappa_ff using the macro treatment. For simple mode, we
   * can simply use radiation which **SHOULD** return the continuum opacity
   * as well, plus something from induced Compton heating. In either cases,
   * we still then need to add the optical depth from electron scattering at
   * the end.
   */

  kappa_tot = 0;

  if (geo.rt_mode == RT_MODE_2LEVEL)
  {
    kappa_tot += radiation (pextract, smax);
  }
  else if (geo.rt_mode == RT_MODE_MACRO)
  {
    if (wind_cell->vol > 0)
    {
      kappa_tot += kappa_bf (plasma_cell, freq_inner, 0);
      kappa_tot += kappa_ff (plasma_cell, freq_inner);
    }
  }

  kappa_tot += klein_nishina (mean_freq) * plasma_cell->ne * zdom[ndom].fill;

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
tau_extract (WindPtr w, PhotPtr porig, double *col_den, double *tau)
{
  int ierr;
  int ndom;
  int istat;
  int n_trans_space;
  int const max_trans_space = 10;
  double norm[3];
  struct photon pextract;

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
        Error ("%s:%s:%i: tau_extract photon transport ended due to too many translate_in_space\n", __FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
      }
    }
    else if ((pextract.grid = where_in_grid (ndom, pextract.x)) >= 0)   // Move the photon in the wind
    {
      ierr = calculate_tau (w, &pextract, col_den, tau);
      if (ierr)
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
    Error ("%s:%s:%i: photon hit central source or disk incorrectly istat = %i\n", __FILE__, __func__, __LINE__, istat);
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
  int wind_index, plasma_index;
  int static xloc_init = 0;

  double x_loc;
  double static upward_x_loc;
  double plasma_density, max_plasma_density;

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

    xloc_init++;
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
 * the surface of the central source, taking into account the current direction
 * of the sight-line being extracted. It's currently not possible, without a
 * tedious re-write, to place the photon at the very origin of the grid when
 * there is a central source because of how we check boundries.
 *
 * Note that photons are initialised with a weight of f_tot as photons are
 * required to have weight, but since functions do not care about the weight of
 * the photon, it is set to something large to make sure it does not get
 * destroyed by accident somehow.
 *
 * ************************************************************************** */

int
create_tau_diag_phot (PhotPtr pout, double nu, double *lmn)
{
  double theta;

  if (nu < 0)
  {
    Error ("%s:%s:%i: photon can't be created with negative nu\n", __FILE__, __func__, __LINE__);
    return EXIT_FAILURE;
  }

  stuff_v (lmn, pout->lmn);

  pout->freq = pout->freq_orig = nu;
  pout->origin = pout->origin_orig = PTYPE_DISK;
  pout->istat = P_INWIND;
  pout->w = pout->w_orig = geo.f_tot;
  pout->tau = 0.0;

  theta = acos (lmn[2]);
  pout->x[0] = geo.rstar + EPSILON;
  pout->x[1] = 0.0;
  pout->x[2] = geo.rstar * sin (theta) + EPSILON;

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
  double default_angles[] = { 0.0, 10.0, 30.0, 45.0, 60.0, 75.0, 85.0, 90.0 };
  int const n_default_angles = sizeof default_angles / sizeof default_angles[0];

  /*
   * Use the angles specified for by the user for spectrum generation, this
   * requires for xxspec to be initialised.
   */

  if (geo.nangles && xxspec != NULL)
  {
    Log_silent ("Using inclinations provided for spectrum cycles as defined in xxspec\n");

    N_ANGLES = geo.nangles;
    OBSERVERS = calloc (geo.nangles, sizeof *OBSERVERS);

    if (OBSERVERS == NULL)
    {
      memory_req = geo.nangles * sizeof *OBSERVERS;
      Error ("%s:%s:%i: cannot allocate %d bytes for observers array\n", __FILE__, __func__, __LINE__, memory_req);
      Exit (1);
    }
    else
    {
      for (iangle = MSPEC; iangle < MSPEC + geo.nangles; iangle++)
      {
        strcpy (OBSERVERS[iangle - MSPEC].name, xxspec[iangle].name);
        stuff_v (xxspec[iangle].lmn, OBSERVERS[iangle - MSPEC].lmn);
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
    OBSERVERS = calloc (n_default_angles, sizeof *OBSERVERS);

    if (OBSERVERS == NULL)
    {
      memory_req = n_default_angles * sizeof *OBSERVERS;
      Error ("%s:%s:%i: cannot allocate %d bytes for observers array\n", __FILE__, __func__, __LINE__, memory_req);
      Exit (1);
    }
    else
    {
      for (iangle = 0; iangle < n_default_angles; iangle++)
      {
        sprintf (OBSERVERS[iangle].name, "A%02.0fP%04.2f", default_angles[iangle], default_phase);
        OBSERVERS[iangle].lmn[0] = sin (default_angles[iangle] / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
        OBSERVERS[iangle].lmn[1] = sin (default_angles[iangle] / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
        OBSERVERS[iangle].lmn[2] = cos (default_angles[iangle] / RADIAN);
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
  int ierr;
  int ispec, ifreq;
  int memory_req;

  double tau, column;
  double freq;
  double freq_min, freq_max, dfreq;
  double wave_min, wave_max;

  double *current_observer;
  double *tau_spectrum;

  struct photon ptau;

  Log ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (N_ANGLES * NWAVE, sizeof *tau_spectrum);
  if (tau_spectrum == NULL)
  {
    memory_req = N_ANGLES * NWAVE * sizeof *tau_spectrum;
    Error ("%s:%s:%i: cannot allocate %d bytes for tau_spectrum\n", __FILE__, __func__, __LINE__, memory_req);
    Exit (1);
  }

  /*
   * Define the limits of the spectra in both wavelength and frequency. The
   * size of the frequency and wavelength bins is also defined. Note that the
   * wavelength limits MUST be in units of Angstroms.
   * TODO: make wavelength limits an advanced diag parameter
   */

  wave_min = 100;
  wave_max = 10000;

  freq_max = VLIGHT / (wave_min * ANGSTROM);
  freq_min = VLIGHT / (wave_max * ANGSTROM);
  dfreq = (freq_max - freq_min) / NWAVE;
  kbf_need (freq_min, freq_max);

  Log_silent ("Minimum wavelength : %.0f Angstroms\n", wave_min);
  Log_silent ("Maximum wavelength : %.0f Angstroms\n", wave_max);
  Log_silent ("Wavelength bins    : %i\n\n", NWAVE);

  /*
   * Now create the optical depth spectra for each observer angle
   */

  for (ispec = 0; ispec < N_ANGLES; ispec++)
  {
    Log ("  - Creating spectrum: %s\n", OBSERVERS[ispec].name);

    freq = freq_min;
    current_observer = OBSERVERS[ispec].lmn;

    for (ifreq = 0; ifreq < NWAVE; ifreq++)
    {
      tau = 0;
      column = 0;

      ierr = create_tau_diag_phot (&ptau, freq, current_observer);
      if (ierr == EXIT_FAILURE)
      {
        Log ("%s:%s:%i: skipping photon of frequency %e\n", __FILE__, __func__, __LINE__, freq);
        continue;
      }

      reposition_tau_photon (&ptau);

      ierr = tau_extract (w, &ptau, &column, &tau);
      if (ierr == EXIT_FAILURE)
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
  int itau, ispec;
  int ierr;
  int memory_req;

  double nu;
  double tau, column;
  double *current_observer;
  double *tau_store, *column_store;

  struct photon ptau;

  tau_store = calloc (N_ANGLES * N_TAU, sizeof *tau_store);
  if (tau_store == NULL)
  {
    memory_req = N_ANGLES * N_TAU * sizeof *tau_store;
    Error ("%s:%s:%i: cannot allocate %d bytes for tau_store\n", __FILE__, __func__, __LINE__, memory_req);
    Exit (1);
  }

  column_store = calloc (N_ANGLES, sizeof *tau_store);
  if (column_store == NULL)
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
    current_observer = OBSERVERS[ispec].lmn;

    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = 0.0;
      column = 0.0;
      nu = TAU_DIAG_OPACS[itau].nu;

      /*
       * Create the tau diag photon and point it towards the extract angle
       */

      ierr = create_tau_diag_phot (&ptau, nu, current_observer);
      if (ierr == EXIT_FAILURE)
      {
        Log ("%s:%s:%i: skipping photon of frequency %e\n", __FILE__, __func__, __LINE__, nu);
        tau_store[ispec * N_TAU + itau] = -1;
        column_store[ispec] = -1;
        continue;
      }

      reposition_tau_photon (&ptau);

      ierr = tau_extract (w, &ptau, &column, &tau);
      if (ierr == EXIT_FAILURE)
      {
        column = -1;
        tau = -1;
      }

      tau_store[ispec * N_TAU + itau] = tau;
      column_store[ispec] = column;
    }
  }

  print_tau_angles (tau_store, column_store);

  free (tau_store);
  free (column_store);
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
   * We have to set this flag, so radiation returns without updating the
   * radiation field and etc estimators. Don't worry, we switch it back to
   * what it was later...
   */

  geo.ioniz_or_extract = 0;

  init_tau_diag_angles ();
  tau_integrate_angles (w);
  create_tau_spectrum (w);

  if (OBSERVERS != NULL)
    free (OBSERVERS);

#ifdef MPI_ON
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  geo.ioniz_or_extract = 1;

  xsignal (files.root, "%-20s Optical depth diagnostics completed\n", "TAU");
  Log (" Completed optical depth diagnostics. The elapsed TIME was %f\n", timer ());
}

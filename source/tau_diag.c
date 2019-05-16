/** ************************************************************************* */
/**
 * @file     tau_diag.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    The main resting place for all functions related to providing
 *           optical depth diagnostic.
 *
 * @details
 *
 * HOW TO ADD MORE OPTICAL DEPTHS
 * ------------------------------
 *
 * In theory it should be relatively straight forwards to add more optical
 * depths if desired. One should only need to update the the two arrays named
 * TAU_NAME and PHOT_FREQ.
 *
 * All one should be required to do is to add a name for the optical depth in
 * TAU_NAME and to append the relevant photon frequency in PHOT_FREQ. This is
 * done as the function radiation is called to calculate the opacity for a
 * photon of a given frequency in the current cell. TAU_NAME is merely just a
 * way to give a name to the optical depth for display purposes.
 *
 * If one wishes to add a "special" opacity, which is an opacity which will
 * require it's own function to calculate, such as the Rosseland mean opacity,
 * then one should update the enumerator SPEC_OPAC with a new label for that
 * opacity, as well as the PHOT_FREQ and TAU_NAME arrays as described above.
 * However, PHOT_FREQ should be given a distinctive value as in the function
 * tau_diag (), the value of opac_type will need to be set as this will be
 * required in the function find_tau () to calculate the correct opacity. Thus,
 * one  should also add another else if branch for the variable opac_type in
 * find_tau ().
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

#include <stdio.h>
#include <string.h>

#include "atomic.h"
#include "python.h"

/*
 * SPEC_OPAC:
 * Enumerator to label any special opacities which are used. This is mainly for
 * the use of the Rosseland and Planck mean opacities.
 */

enum SPEC_OPAC
{
  ROSSEL_MEAN = 1,
  PLANCK_MEAN = 2,
};

/*
 * TAU_NAME:
 * This array is used for book keeping purposes when printing information to
 * the screen.
 * NOTE: this array is global for convenience purposes if more optical depths
 *       wish to be added
 */

char *TAU_NAME[] = {
  "Rosseland",
  "Planck",
  "Lyma_850A",
  "Bal_3600A",
};

/*
 * PHOT_FREQ:
 * This array is used to set the frequency for a tau diag photon. By setting the
 * photon's frequency, Python should be able to calculate the frequency
 * dependent opacity for that photon.
 * NOTE: this array is global for convenience purposes if more optical depths
 *       wish to be added
 */

double PHOT_FREQ[] = {
  -ROSSEL_MEAN,
  -PLANCK_MEAN,
  3.526970e+15,
  8.327568e+14,
};

#define N_TAU (sizeof (PHOT_FREQ) / sizeof (*PHOT_FREQ))
#define N_ANGLES (nspectra - MSPEC)
#define MAX_TRANS_SPACE 10

/* ************************************************************************* */
/**
 * @brief           Print the various optical depths calculated using this
 *                  routine
 *
 * @param[in]       double   tau_store    The 2d array containing the optical
 *                                        depth for each observer and tau
 *
 * @return          void
 *
 * @details
 *
 * Simply prints the different optical depths for each angle and optical depth.
 *
 * ************************************************************************** */

void
print_tau_table (double tau_store[N_ANGLES][N_TAU])
{
  int itau;
  int ispec;
  int line_len;
  int const MAX_COL = 120;

  double tau;

  char tmp_str[50];
  char observer_name[40];

  Log ("\nOptical depths along the defined line of sights:\n\n");

  for (ispec = MSPEC; ispec < nspectra; ispec++)
  {
    strcpy (observer_name, xxspec[ispec].name);
    Log ("%s\n--------\n", observer_name);

    line_len = 0;
    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = tau_store[ispec - MSPEC][itau];
      line_len += sprintf (tmp_str, "tau_%-9s: %3.2e  ", TAU_NAME[itau], tau);

      if (line_len > MAX_COL)
      {
        line_len = 0;
        Log ("\n\t");
      }

      Log ("%s", tmp_str);
    }

    Log ("\n\n");
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
 * @param[in,out]   double    *tau        The optical depth across the cell
 *
 * @return          int       istat       The current photon status
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
 * TODO: implement Rosseland mean opacity
 * TODO: implement Planck mean opacity
 *
 * ************************************************************************** */

int
find_tau (WindPtr w, PhotPtr pextract, int opac_type, double *tau)
{
  int istat;

  double smax;
  double kappa_tot;

  if ((pextract->grid = where_in_grid (w[pextract->grid].ndom, pextract->x)) < 0)
  {
    Error ("%s:%i:%s: pextract is not in grid\n", __FILE__, __LINE__, __func__);
    return pextract->grid;
  }

  if ((smax = find_smax (pextract)) < 0)
  {
    Error ("%s:%i:%s: abnormal value of smax for photon\n", __FILE__, __LINE__, __func__);
    return -1;
  }

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
 * @param[out]      double    *tau       The optical depth from photon origin to
 *                                       the observer
 *
 * @return          void
 *
 * @details
 *
 * The direction of the observer is set in porig.lmn and should be set prior
 * to passing a photon to this function. If any values of lmn are negative, the
 * phton will translate in the negative direction. However, have no fear as this
 * is normal and is fine due to the assumed symmetry of models in Python.
 *
 * ************************************************************************** */

void
extract_tau (WindPtr w, PhotPtr porig, int opac_type, double *tau)
{
  int ndom;
  int istat;
  int n_trans_space;

  double norm[3];

  struct photon pextract;

  // TODO: remove when Planck and Rosseland mean have been implemeneted properly
  if (opac_type == ROSSEL_MEAN || opac_type == PLANCK_MEAN)
    return;

  istat = P_INWIND;             // assume photon is in wind for initialisation reasons
  stuff_phot (porig, &pextract);

  /*
   * In the case of an observer which does not look into the wind from the
   * photon origin, the while loop would become infinite as translate_in_space
   * would fail to find the wind hence a limit is imposed here for the number
   * of attempts to find the nearest wind boundary
   */

  n_trans_space = 0;

  while (istat == P_INWIND)
  {
    /*
     * The first part of this while loop is concerned with translating the
     * photon either in space, i.e. in the non-wind grid cells, or translates
     * the photon one grid cell at a time using find_tau.
     */

    if (where_in_wind (pextract.x, &ndom) < 0)
    {
      translate_in_space (&pextract);
      if (++n_trans_space > MAX_TRANS_SPACE)
      {
        Error ("%s:%i:%s: photon transport ended due to too many translate_in_space\n", __FILE__, __LINE__, __func__);
        break;
      }
    }
    else if ((pextract.grid = where_in_grid (ndom, pextract.x)) >= 0)
    {
      istat = find_tau (w, &pextract, opac_type, tau);
      if (istat != P_INWIND)
        break;
    }
    else
    {
      pextract.istat = -1;
      Error ("%s:%i:%s: photon in unknown location,  grid stat %i\n", __FILE__, __LINE__, __func__, pextract.grid);
      return;
    }

    /*
     * Now that the photon has been translated, we updated istat using the walls
     * function which should return an istat which is not P_INWIND when the
     * photon has exited in the wind -- hopefully when it has escaped
     */

    istat = walls (&pextract, porig, norm);
  }

  if (istat == P_HIT_STAR || istat == P_HIT_DISK)
    Error ("%s:%i:%s: photon hit star or disk, istat = %i\n", __FILE__, __LINE__, __func__, istat);

  if (walls (&pextract, porig, norm) == -1)
    Error ("%s:%i:%s: abnormal return from walls for photon\n", __FILE__, __LINE__, __func__);
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
 * @return          void
 *
 * @details
 *
 * This routine assumes that the user wishes for the photon to be generated from
 * the inner-most radius of the disk or from the central source of the
 * simulation. Note that photons are initialised with a weight of f_tot as
 * photons are required to have weight, but since the diagnostic functions do
 * not care about the weight of the photon, it is set to something large.
 *
 * TODO: add more flexible photon placement: placement of the photon is v v important
 *
 * ************************************************************************** */

void
tau_diag_phot (PhotPtr pout, double nu)
{
  pout->freq = pout->freq_orig = nu;
  pout->origin = pout->origin_orig = PTYPE_DISK;
  pout->istat = P_INWIND;
  pout->w = pout->w_orig = geo.f_tot;
  pout->x[0] = geo.rstar;
  pout->x[1] = 0;
  pout->x[2] = geo.rstar;
  pout->tau = 0;
}

/* ************************************************************************* */
/**
 * @brief     The main controlling function for the optical depth diagnostics.
 *
 * @return    void
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
 *  - Average freq for the photon bands
 *
 * The aim of these diagnostic numbers it to provide some sort of metric on the
 * optical thickness of the current model.
 *
 * TODO: fixed angles rather than angles in xxspec could be better - but this is complicated for binary systems
 * TODO: method for calling this for individual cells instead of the entire wind
 *
 * ************************************************************************** */

void
tau_diag (WindPtr w)
{
  int itau;
  int ispec;
  int opac_type;

  double nu;
  double tau;
  double *observer;

  struct photon ptau;

  /*
   * EP:
   * 2d array for storing the optical depths for each observer angle and tau..
   * I tried assigning this dynamically, in case this gets allocated on stack
   * and causes a segfault, but this for some reason caused segfaults in other
   * parts of the code.
   */

  double tau_store[N_ANGLES][N_TAU];

  /*
   * Switch to extract mode - this is bad generally not great practise but it
   * enables radiation (p, ds) to return the opacity in a cell without updating
   * the monte carlo estimators
   */

  geo.ioniz_or_extract = 0;

  /*
   * MSPEC here is used as this is where the actual observer angles start
   */

  for (ispec = MSPEC; ispec < nspectra; ispec++)
  {
    observer = xxspec[ispec].lmn;
    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = 0;
      nu = PHOT_FREQ[itau];

      /*
       * If the frequency is negative, a mean opacity will be used instead
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
        opac_type = N_TAU;

      /*
       * Create the tau diag photon and point it towards the extract angle
       */

      tau_diag_phot (&ptau, nu);
      stuff_v (observer, ptau.lmn);

    /* **********************************************************************

          It is a great issue with me on where to place the extracted photon.
          This bit of code was an experiment to see if I should place the
          photon at -Rstar if lmn[0] was zero.. but the results were very
          sensitive to the placement. Still unsure which is better right now.

          if (ptau.lmn[0] < 0)
            ptau.x[0] *= -1;

       ********************************************************************** */

      /*
       * Extract tau and add it to the tau storage array
       */

      extract_tau (w, &ptau, opac_type, &tau);
      tau_store[ispec - MSPEC][itau] = tau;
    }
  }

  print_tau_table (tau_store);

  /*
   * Switch back to ionization mode - may be redundant but if this is called
   * during each ionisation cycle then this will be required...
   */

  geo.ioniz_or_extract = 1;
}

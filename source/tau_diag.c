/** ************************************************************************* */
/**
 * @file     tau_diag.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    The main resting place for all functions related to providing
 *           optical depth diagnostic.
 *
 * ************************************************************************** */

#include <stdio.h>
#include <string.h>

#include "atomic.h"
#include "python.h"


#define MAX_COL 120
#define MAX_TRANS_SPACE 10
#define N_ANGLES nspectra - MSPEC

enum OPACITIES                  // TODO: likely redundant, could use sizeof to automatically find N_TAU
{
  ROSSELAND,
  PLANCK,
  LYMAN_EDGE,                   // 912 A
  BALMER_EDGE,                  // 3646 A
  N_TAU                         // Used as a counter for the number of taus available
};

/*
 * This array is for book keeping purposes when printing information to the
 * screen.
 */

char *OPACITY_NAMES[] = {
  "Rosseland",
  "Planck",
  "Lyman",
  "Balmer"
};

/*
 * This array is used to set the frequency a tau diag photon should be so the
 * opacity in each cell will be calculated for a photon of this frequency.
 */

double PHOTON_FREQS[] = {
  -1                            // TODO: not sure what to set mean opacity frequencies to...
    - 2,
  3.28719e15,
  8.22503e14,
};


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
 * TODO: modularise in future so w could also be a single cell
 * TODO: implement Rosseland mean
 * TODO: implement Planck mean
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
    Error ("find_tau (%s:%i): pextract is not in grid\n", __FILE__, __LINE__);
    return pextract->grid;
  }

  /*
   * smax is the distance a photon can move in the cell in its current direction
   * until it hits the boundary of the cell
   */

  smax = find_smax (pextract);
  if (smax < 0)
    return -1;

  if (opac_type == ROSSELAND)   // Not implemented
  {
    move_phot (pextract, geo.rmax);
    return P_ESCAPE;
  }
  else if (opac_type == PLANCK) // Not implemented
  {
    move_phot (pextract, geo.rmax);
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
 * @param[in]       WindPtr   w          A pointer to the wind
 * @param[in]       PhotPtr   porig      The photon packet to extract
 * @param[in]       int       opac_type  An indicator of whether to use a mean
 *                                       opacity or not.
 * @param[out]      double    *tau       The optical depth from photon origin to
 *                                       the observer
 *
 * @return          void
 *
 * @details
 *
 * The direction of the observer is set in porig.lmn and should be set prior
 * to passing a photon to this function. If any values of lmn are negative, the
 * photon will translate in the negative direction. However, have no fear as this
 * is normal and is fine due to assumed symmetries.
 *
 * TODO: modularise in future so w could also be a single cell
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

  istat = P_INWIND;
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
        Error ("extract_tau (%s:%i): photon transport ended due to too many translate_in_space\n", __FILE__, __LINE__);
        return;
      }
    }
    else if ((pextract.grid = where_in_grid (ndom, pextract.x)) >= 0)
    {
      pextract.istat = istat = find_tau (w, &pextract, opac_type, tau);
      if (istat == -1)
      {
        Error ("extract_tau (%s:%i): abnormal value of smax found in find_tau\n", __FILE__, __LINE__);
        return;
      }
    }
    else
    {
      pextract.istat = -1;
      Error ("extract_tau (%s:%i): photon in unknown location,  grid stat %i\n", __FILE__, __LINE__, pextract.grid);
      return;
    }

    /*
     * Now that the photon has been translated, we updated istat using the walls
     * function which should return an istat which is not P_INWIND when the
     * photon has exited in the wind -- hopefully when it has escaped
     */

    istat = walls (&pextract, porig, norm);
  }

  if (walls (&pextract, porig, norm) == -1)
    Error ("extract_tau (%s: %i): abnormal return from walls for photon\n", __FILE__, __LINE__);
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
 * TODO: add more flexible photon placement
 *
 * ************************************************************************** */

void
create_tau_diag_phot (PhotPtr pout, double nu)
{
  pout->freq = pout->freq_orig = nu;
  pout->origin = pout->origin_orig = PTYPE_DISK;
  pout->istat = P_INWIND;
  pout->w = pout->w_orig = geo.f_tot;
  pout->x[0] = geo.rstar;
  pout->x[1] = 0.0;
  pout->x[2] = 0.0;
  pout->tau = 0;
}

/* ************************************************************************* */
/**
 * @brief           Print the various optical depths calculated using this
 *                  routine
 *
 * @param[in]       double   **tau_store  The 2d array containing the optical
 *                                        depth for each observer and tau
 *
 * @return          void
 *
 * @details
 *
 * ************************************************************************** */

void
print_tau_table (double tau_store[N_ANGLES][N_TAU])
{
  int itau;
  int ispec;
  int line_len;

  double tau;

  char tmp_str[50];
  char observer_name[40];

  Log ("\nOptical depths along the observer's line of sight:\n\n");

  for (ispec = MSPEC; ispec < nspectra; ispec++)
  {
    strcpy (observer_name, xxspec[ispec].name);
    Log ("%s\n\t", observer_name);

    line_len = 0;
    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = tau_store[ispec - MSPEC][itau];
      line_len += sprintf (tmp_str, "tau_%-9s: %3.2e  ", OPACITY_NAMES[itau], tau);
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
 * @brief           The main steering function for the optical depth
 *                  diagnostics.
 *
 * @param[in]       WindPtr     A pointer to the wind
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
 *  - Average freq for the photon bands
 *
 * The aim of these diagnostic numbers it to provide some sort of metric on the
 * optical thickness of the current model.
 *
 * TODO: fixed angles rather than the angles in xxspec would be better
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
   * I wanted to change this to dynamic allocation or to a 1d array to make it
   * less "ugly", but I came across issues with calloc memory corruptions and
   * weird segfaults elsewhere in the code (related to OOE?).
   */

  double tau_store[N_ANGLES][N_TAU];

  /*
   * Switch to extract mode - this is bad generally not great practise but it
   * enables radiation (p, ds) to return the opacity in a cell without updating
   * the monte carlo estimators
   */

  geo.ioniz_or_extract = 0;

  // MSPEC here is used as this is where the actual observer angles start
  for (ispec = MSPEC; ispec < nspectra; ispec++)
  {
    observer = xxspec[ispec].lmn;
    for (itau = 0; itau < N_TAU; itau++)
    {
      tau = 0;
      nu = PHOTON_FREQS[itau];

      // If the frequency is negative, a mean opacity will be used instead :-)
      if (nu == -1)
      {
        opac_type = ROSSELAND;
        nu = 1e15;
      }
      else if (nu == -2)
      {
        opac_type = PLANCK;
        nu = 1e15;
      }
      else
        opac_type = N_TAU;

      // Create the tau diag photon and point it towards the extract angle
      create_tau_diag_phot (&ptau, nu);
      stuff_v (observer, ptau.lmn);

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

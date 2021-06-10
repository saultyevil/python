
/***********************************************************/
/** @file  gradv.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Routines related to calculating velocity gradients
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief      find the gradient velocity vector v, dv_ds for a photon
 * at a certain position travelling in a certain direction in the cmf
 * frame.
 *
 * @param [in] PhotPtr  p   A photon
 * @return     dvds    on successful completion
 * the value of where in grid if the photon outside the wind grid.
 *
 * @details
 *
 * ### Notes ###
 * For spherical coordinates the routine calculates the gradient
 * on the fly, that is moves the photon a bit and calculates dvds
 *
 * For 2d systems, the velocity gradient is calculated using
 * the velocity gradient tensors, which contain the velocity
 * gradient in the xz plane.
 *
 * It's not clear how much faster the method used for 2d systems
 * actually is, and one might want to consider the on-the-
 * fly appoach for all systems.
 *
 **********************************************************/

double
dvwind_ds_cmf (p)
     PhotPtr p;
{
  double v_grad[3][3];
  double lmn[3], dvel_ds[3], dvds;
  int j, k, nn;
  struct photon pp;
  int nnn[4], nelem;            // At present the largest number of dimensions in the grid is 2
  double frac[4];
  double x;
  int ndom;

  ndom = wmain[p->grid].ndom;

  /* We want the change in velocity along the line of sight, but we
     need to be careful because of the fact that we have elected to
     combine the upper and lower hemispheres in the wind array.  Since
     we are only concerned with the the scalar dv_ds, the safest thing
     to do is to create a new photon that is only in the upper hemisphere
     02jan ksl */

  stuff_phot (p, &pp);
  if (pp.x[2] < 0.0)
  {                             /*move the photon to the northern hemisphere */
    pp.x[2] = -pp.x[2];
    pp.lmn[2] = -pp.lmn[2];
  }

  /* JM 1411 -- ideally, we want to do an interpolation on v_grad here. However,
     the interpolation was incorrect in spherical coordinates (see issue #118).
     For the moment, I've adopted an on the fly method for spherical coordinates.
     This should be improved by figuring out out how the velocity gradient
     tensor ought to be rotated in order to give the right answer for spherical
     coordinates. */

  if (zdom[ndom].coord_type == SPHERICAL || USE_GRADIENTS == FALSE)
  {
    struct photon pnew;
    double v1[3], v2[3], dv[3], diff[3];
    double ds;
    /* choose a small distance which is dependent on the cell size */
    vsub (pp.x, wmain[pp.grid].x, diff);
    ds = 0.001 * length (diff);
    /* calculate the velocity at the position of the photon */
    /* note we use model velocity, which could potentially be slow,
       but avoids interpolating (see #118) */
    model_velocity (ndom, pp.x, v1);

    /* copy the photon and move it by ds, and evaluate the velocity
       at the new point */
    stuff_phot (&pp, &pnew);
    /*
     * Put the photon into the observer frame, this way move_phot won't throw
     * an error. Should be ok since we are only using move_photon to move a
     * photon some vector ds.
     * */
    pnew.frame = F_OBSERVER;
    move_phot (&pnew, ds);
    model_velocity (ndom, pnew.x, v2);

    /* calculate the relevant gradient */
    if (rel_mode == REL_MODE_FULL)
    {
      observer_to_local_frame_velocity (v2, v1, dv);
      dvds = length (dv) / ds;
    }
    else
    {
      dvds = fabs (dot (v1, pp.lmn) - dot (v2, pp.lmn)) / ds;
    }
  }

  else                          // for non spherical coords we interpolate on v_grad
  {

    coord_fraction (ndom, 0, pp.x, nnn, frac, &nelem);


    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        x = 0;
        for (nn = 0; nn < nelem; nn++)
          x += wmain[nnn[nn]].v_grad[j][k] * frac[nn];

        v_grad[j][k] = x;

      }
    }

    /* v_grad is in cylindrical coordinates, or more precisely intended
       to be azimuthally symmetric.  One could either
       (a) rotate  v_grad to be correct at the position of the photon or
       (b) rotate the direction of photon travel so that is is correct
       (assuming azimuthal symmetry) in the xz plane.

       Possibility b is more straightforward and that is what is done
     */

    project_from_xyz_cyl (pp.x, pp.lmn, lmn);

    dvds = dot_tensor_vec (v_grad, lmn, dvel_ds);

    /* Note that the vector dvel_ds is also in an azimuthally symmetric system in the
     * xx plane, and could be rotated back if it were needed
     */
  }

  if (sane_check (dvds))
  {
    Error ("dvwind_ds_cmf: sane_check %f\n", dvds);
  }

  return (dvds);

}



#define N_DVDS_AVE	10000

/**********************************************************/
/**
 * @brief      Calculate the direction averaged (and maximum)
 * dv_ds in each grid cell of the wind
 *
 * @return     Always returns 0
 *
 * @details
 * The routine cycles through all of the cells in the wind, and calculates
 * the aveage value of dv_ds at the center of the wind cell by randomly
 * generating directions and then calculating dv_ds in these directions
 *
 *
 *
 * ### Notes ###
 * The routine is called during the intialization process and fills
 * the following elements of wmain
 *
 *  * dvds_ave - the average dvds
 *
 * There is an advanced mode which prints this information to
 * file.
 *
 * 210303 - ksl - Removed calculation of dvds_max from this 
 * routine because we want this at the edges of cells, so one 
 * can interpolate.  Portions of the old routine still remain.
 **********************************************************/


int
dvds_ave ()
{
  struct photon p, pp;
  double v_zero[3], delta[3], vdelta[3], diff[3];
  double sum, dvds, ds;
  double dvds_max, lmn[3];
  int n;
  int icell;
  double dvds_min, lmn_min[3];
  char filename[LINELENGTH];
  int ndom;


  /* Open a diagnostic file if print_dvds_info is non-zero */
  strcpy (filename, files.root);
  strcat (filename, ".dvds.diag");
  if (modes.print_dvds_info)
  {
    optr = fopen (filename, "w");
  }

  for (icell = 0; icell < NDIM2; icell++)
  {
    ndom = wmain[icell].ndom;

    dvds_max = 0.0;             // Set dvds_max to zero for the cell.
    dvds_min = 1.e30;           // TEST


    /* Find the center of the cell */

    stuff_v (wmain[icell].xcen, p.x);

    /* Define a small length */

    vsub (p.x, wmain[icell].x, diff);
    ds = 0.001 * length (diff);

    /* Find the velocity at the center of the cell */
    model_velocity (ndom, p.x, v_zero);

    sum = 0.0;
    for (n = 0; n < N_DVDS_AVE; n++)
    {
      randvec (delta, ds);
      if (p.x[2] + delta[2] < 0)
      {                         // Then the new position would punch through the disk
        delta[0] = (-delta[0]); // So we reverse the direction of the vector
        delta[1] = (-delta[1]);
        delta[2] = (-delta[2]);
      }
      vadd (p.x, delta, pp.x);
      model_velocity (ndom, pp.x, vdelta);
      vsub (vdelta, v_zero, diff);
      dvds = length (diff);

      /* Find the maximum and minimum values of dvds and the direction
       * for this
       */

      if (dvds > dvds_max)
      {
        dvds_max = dvds;
        renorm (delta, 1.0);
        stuff_v (delta, lmn);
      }
      if (dvds < dvds_min)
      {
        dvds_min = dvds;
        renorm (delta, 1.0);
        stuff_v (delta, lmn_min);
      }

      sum += dvds;

    }

    /* Store the results in wmain */
    wmain[icell].dvds_ave = sum / (N_DVDS_AVE * ds);

    if (modes.print_dvds_info)
    {
      fprintf (optr,
               "%d %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e \n",
               icell, p.x[0], p.x[1], p.x[2], dvds_max / ds,
               dvds_min / ds, lmn[0], lmn[1], lmn[2], lmn_min[0], lmn_min[1], lmn_min[2], dot (lmn, lmn_min));
    }

  }


  if (modes.print_dvds_info)
    fclose (optr);

  return (0);
}



/**********************************************************/
/**
 * @brief      Calculate the maximum
 * dv_ds in each grid cell of the wind
 *
 * @return     Always returns 0
 *
 * @details
 * The routine cycles through all of the cells in the wind, and calculates
 * the maxium value of dv_ds at the corner of the wind cell by randomly
 * generating directions and then calculating dv_ds in these directions
 *
 *
 * ### Notes ###
 * The routine is called during the intialization process and fills
 * the following elements of wmain
 *
 *  * dvds_max - the maximum value of dvds
 *
 * Unlike dvds_ave, these values are at the corners of cells,
 * and are intended to be interpolated.
 *
 * There is an advanced mode which prints this information to
 * file.
 **********************************************************/


int
dvds_max ()
{
  struct photon p, pp;
  double v_zero[3], delta[3], vdelta[3], diff[3];
  double sum, dvds, ds;
  double dvds_max, lmn[3];
  int n;
  int icell;
  double dvds_min, lmn_min[3];
  char filename[LINELENGTH];
  int ndom;


  /* Open a diagnostic file if print_dvds_info is non-zero */
  strcpy (filename, files.root);
  strcat (filename, ".dvds.diag");
  if (modes.print_dvds_info)
  {
    optr = fopen (filename, "w");
  }

  for (icell = 0; icell < NDIM2; icell++)
  {
    ndom = wmain[icell].ndom;

    dvds_max = 0.0;
    dvds_min = 1.e30;


    /*  x is the corner of the cell */

    stuff_v (wmain[icell].x, p.x);

    /* Define a small length */

    vsub (wmain[icell].xcen, p.x, diff);
    ds = 0.000001 * length (diff);

    /* Find the velocity at the center of the cell */
    model_velocity (ndom, p.x, v_zero);

    sum = 0.0;
    for (n = 0; n < N_DVDS_AVE; n++)
    {
      randvec (delta, ds);
      if (p.x[2] + delta[2] < 0)
      {                         // Then the new position would punch through the disk
        delta[0] = (-delta[0]); // So we reverse the direction of the vector
        delta[1] = (-delta[1]);
        delta[2] = (-delta[2]);
      }
      vadd (p.x, delta, pp.x);
      model_velocity (ndom, pp.x, vdelta);
      vsub (vdelta, v_zero, diff);
      dvds = length (diff);

      /* Find the maximum and minimum values of dvds and the direction
       * for this
       */

      if (dvds > dvds_max)
      {
        dvds_max = dvds;
        renorm (delta, 1.0);
        stuff_v (delta, lmn);
      }
      if (dvds < dvds_min)
      {
        dvds_min = dvds;
        renorm (delta, 1.0);
        stuff_v (delta, lmn_min);
      }


    }

    /* Store the results in wmain */
    wmain[icell].dvds_max = dvds_max / ds;

    if (modes.print_dvds_info)
    {
      fprintf (optr,
               "%d %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e \n",
               icell, p.x[0], p.x[1], p.x[2], dvds_max / ds,
               dvds_min / ds, lmn[0], lmn[1], lmn[2], lmn_min[0], lmn_min[1], lmn_min[2], dot (lmn, lmn_min));
    }

  }


  if (modes.print_dvds_info)
    fclose (optr);

  return (0);
}




/**********************************************************/
/**
 * @brief      Calculate the maximum
 * dv_ds at a particular position in a grid
 *
 * @param [in] PhotPtr  p   A photon

 * @return     Returns dvds_max at the position of the photon
 *
 * @details
 * The routine interpolates dvds_max given the position of
 * a photon in a cell
 *
 * dvds_max at the vertex points of cells must have been
 * initialized using the routine dvds_max
 *
 * ### Notes ###
 *
 * The routine uses both the position of the photon and 
 * the grid cell in which the photon exitsts, so this
 * must be acurrate.
 **********************************************************/


double
get_dvds_max (p)
     PhotPtr p;
{
  int ndom, nn, nnn[4], nelem;
  double frac[4];
  double dvds;

  ndom = wmain[p->grid].ndom;

  coord_fraction (ndom, 0, p->x, nnn, frac, &nelem);

  dvds = 0;

  for (nn = 0; nn < nelem; nn++)
  {
    dvds += wmain[nnn[nn]].dvds_max;
  }

  return dvds;

}

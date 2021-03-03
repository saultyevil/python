
/***********************************************************/
/** @file   setup_disk.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Read parameters that define a disk
 *
 * File containing reverberation mapping functions.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/**********************************************************/
/**
 * @brief       get the parameters needed to define a disk
 *
 * @return Always returns 0
 *
 * Read the parameters, such as the type of disk, the
 * temperature profile, that define a disk
 *
 * The parameters fill variables defined in the geo
 * data structure.
 *
 * ###Notes###
 *
***********************************************************/


int
get_disk_params (void)
{
  char values[LINELENGTH], answer[LINELENGTH];

  rdpar_comment ("Parameters for the Disk (if there is one)");

  if (geo.system_type == SYSTEM_TYPE_STAR)
  {
    strcpy (answer, "none");
  }
  else
  {
    strcpy (answer, "flat");
  }

  sprintf (values, "%d,%d,%d", DISK_NONE, DISK_FLAT, DISK_VERTICALLY_EXTENDED);
  geo.disk_type = rdchoice ("Disk.type(none,flat,vertically.extended)", values, answer);

  if (geo.disk_type == DISK_NONE)
  {
    geo.disk_radiation = 0;
    geo.diskrad = 0;
    return (0);
  }

  strcpy (answer, "yes");
  geo.disk_radiation = rdchoice ("Disk.radiation(yes,no)", "1,0", answer);
  if (geo.disk_radiation)
    get_spectype (geo.disk_radiation, "Disk.rad_type_to_make_wind(bb,models)", &geo.disk_ion_spectype);

  /*
   * Set a default for disk radius for the current system, this is set for a
   * CV or AGN/BH but it set to zero for a star.
   */

  if (geo.system_type == SYSTEM_TYPE_CV)
  {
    geo.diskrad = diskrad (geo.mstar, geo.m_sec, geo.period);
  }
  else if (geo.system_type == SYSTEM_TYPE_AGN || geo.system_type == SYSTEM_TYPE_BH)
  {
    geo.diskrad = 100. * geo.rstar;
  }
  else
  {
    geo.diskrad = 0;
  }

  /*
   * At this point, we query the user on how to set the effective temperature
   * profile for the accretion disc. If the user wants a standard or near-Eddington
   * critical disc, then they are asked for the accretion rate and radius of this
   * disc. Otherwise, they need to provide the filename for the temperature
   * structure and the disc radius will be inferred from this - the accretion
   * rate should therefore not be required.
   */

  geo.disk_tprofile = DISK_TPROFILE_STANDARD;
  strcpy (answer, "standard");
  sprintf (values, "%d,%d,%d", DISK_TPROFILE_STANDARD, DISK_TPROFILE_EDDINGTON_CRITICAL, DISK_TPROFILE_READIN);
  geo.disk_tprofile = rdchoice ("Disk.temperature.profile(standard,eddington,readin)", values, answer);

  if (geo.disk_tprofile == DISK_TPROFILE_STANDARD || geo.disk_tprofile == DISK_TPROFILE_EDDINGTON_CRITICAL)
  {
    /* Query the accretion rate first */

    geo.disk_mdot /= (MSOL / YR);       // Convert to msol/yr to simplify input
    rddoub ("Disk.mdot(msol/yr)", &geo.disk_mdot);
    geo.disk_mdot *= (MSOL / YR);

    /* Now query the radius of the accretion disc */

    rddoub ("Disk.radmax(cm)", &geo.diskrad);
    Log ("geo.diskrad  %e\n", geo.diskrad);
    geo.diskrad_sq = geo.diskrad * geo.diskrad;
  }
  else if (geo.disk_tprofile == DISK_TPROFILE_READIN)
  {
    rdstr ("Disk.T_profile_file", files.tprofile);
    geo.diskrad = read_non_standard_disk_profile (files.tprofile);
    geo.disk_mdot = 0;
  }
  else
  {
    Error ("Setup_disk: This should never occur\n");
    Exit (1);
  }

  /* If diskrad <= geo.rstar set geo.disk_type = DISK_NONE to make any disk transparent anyway. */

  if (geo.diskrad < geo.rstar)
  {
    Log ("Disk radius is less than star radius, so assuming no disk)\n");
    geo.disk_type = DISK_NONE;
  }

  /* Get the additional variables need to describe a vertically extended disk */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    rddoub ("Disk.z0(fractional.height.at.diskrad)", &geo.disk_z0);
    rddoub ("Disk.z1(powerlaw.index)", &geo.disk_z1);
  }

  return (0);
}

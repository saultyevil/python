/* ************************************************************************** */
/**
 * @file  import.h
 * @author EJP
 * @date   Feb 2019
 *
 * @brief    Global structures for containing imported models
 *
 * @details
 * ************************************************************************** */

#define DEFAULT_IMPORT_TEMPERATURE 40000

struct
{
  int ncell;
  int *element;
  double *r;
  double *v_r;
  double *mass_rho;
  double *t_r;

} import_model_1d;

/*
 * i is the row number
 * j is the column number
 */

#define READ_NO_TEMP 9
#define READ_RAD_TEMP 10
#define READ_BOTH_TEMP 11

struct
{
  int ncell;                          // the total number of cells read in
  int ndim, mdim;                     // the number of coordinates in the n and m dimensions
  int *i, *j;                         // the i (row) and j (column) elements
  int *inwind;                        // flag for the cell being inwind or not inwind
  double *x, *z, *r, *theta;          // the x/r or z/theta coordinates of the grid
  double *v_x, *v_y, *v_z;            // the velocity in Cartesian coordinates
  double *mass_rho;                   // the mass density in cgs units
  double *t_e, *t_r;                  // the electron and radiation temperature in Kelvin
  double *wind_x, *wind_z;
  double *wind_midx, *wind_midz;
} import_model_2d;


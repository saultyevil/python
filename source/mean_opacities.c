/** ************************************************************************* */
/**
 * @file     mean_opacities.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    Contains functions used to calculate mean opaicities, such as the
 *           Rosseland mean.
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>

#include "atomic.h"
#include "python.h"

/* ************************************************************************* */
/**
 * @brief     Computes the temperature derivative of the Planck function.
 *
 * @param[in] double nu             The frequency of the black body emission
 * @param[in] double temperature    The temperature of the black body
 *
 * @returns   double dB_dT          The temperature derivative of the Planck
 *                                  function wrt to temperature for a given
 *                                  frequency and temperature
 *
 * @details
 *
 * I take some shortcuts with algebra for dB/dT hence the final expression is a
 * bit odd looking.
 *
 * ************************************************************************** */

double ross_int_temp;

double
dB_nu_dT (double nu)
{
  double x;
  double bnu;
  double dB_dT;

  x = H_OVER_K * nu / ross_int_temp;
  bnu = 2 * H * pow (nu, 3) / pow (C, 2) * 1 / (exp (x) - 1);
  dB_dT = bnu * x / (ross_int_temp * (1 - exp (-x)));

  return dB_dT;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]
 *
 * @return
 *
 * @details
 *
 * ************************************************************************** */

int iter = 0;

double
rosseland_integrand (double nu)
{
  double kappa_es;
  double integrand;

  kappa_es = klein_nishina (nu);
  integrand = 1 / kappa_es * dB_nu_dT (nu);

  return integrand;
}


/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]
 *
 * @return
 *
 * @details
 *
 * ************************************************************************** */

void
set_ross_integ_temp (double temperature)
{
  ross_int_temp = temperature;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]
 *
 * @return
 *
 * @details
 *
 * ************************************************************************** */

void
calculate_rosseland_opacity (void)
{
  int i;
  double integral = 0;
  double ross_const;
  double kappa_ross = 0.0;
  double temperature = 1.5e4;

  double (*f) (double) = &rosseland_integrand;

  Log ("calculating rosseland mean opacity for %e K\n", temperature);

  int n_temps = 1;
  for (i = 0; i < n_temps; i++)
  {
    set_ross_integ_temp (temperature);
    ross_const = PI / (4 * STEFAN_BOLTZMANN * pow (temperature, 3));
    integral = qromb (f, 0, VERY_BIG, EPSILON);
    integral *= ross_const;
    kappa_ross = 1 / integral;
  }

  Log ("Rosseland integral: %e\n", kappa_ross);

}


/* ************************************************************************* */
/**
 * @brief
 *
 * @param[in]
 *
 * @return
 *
 * @details
 *
 * ************************************************************************** */

void
calculate_planck_opacity (void)
{

}

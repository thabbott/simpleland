/**
 * land.c
 *
 * Simple representation of a land surface with prognostic temperature and soil
 * moisture. Documentation for users is in land.h
 *
 * Tristan Abbott // MIT EAPS // September 2018
 */

// Import required interfaces
#include <math.h>
#include <stdlib.h>

// Implements this interface
#include "land.h"

/* Macros */
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))

/* Private constants */
/** Neutral stability Prandtl number */
static const double Pr_0 = 0.74;
/** von Karman constant */
static const double k_vK = 0.40;
/** Bulk Richardson number where unstable solution changes */
static const double Ri_c = -0.2097;
/** Profile function parameters */
static const double alpha_m = 4.7;
static const double beta_m = 15.0;
static const double alpha_h = 6.35135135135;
static const double beta_h = 9.0;

/**
 * Implements land_flux
 */
void simple_land_flux(double rho, double grav,
                      double c_p, double L_v,
                      double theta_s, double theta_atm,
                      double qstar_s, double q_atm,
                      double u_atm, double v_atm,
                      double u_min,
                      double phi, double phi_fc,
                      double phi_pwp, double r_sfc,
		      double z_atm, double z_0,
                      double *shf, double *lhf,
                      double *taux, double *tauy) {

  // Calculate bulk Richardson number
  double u_mag = sqrt(u_atm*u_atm + v_atm*v_atm);
  if (u_mag < u_min) { 
    u_mag = u_min; 
  }
  double Ri_b = grav / theta_s *
    (theta_atm - theta_s) * z_atm / (u_mag * u_mag);

  // Calculate M-O stability parameter
  double zeta;
  double sb;
  double Qb;
  double Pb;
  double thetab;
  double Tb;
  double psi_m;
  double psi_h;
  double x;
  double x_0;
  double y;
  double y_0;
  if (Ri_b >= 0) {
   
    zeta = z_atm / (z_atm - z_0) * log(z_atm / z_0) / 
           (2.0 * alpha_h * (alpha_m * Ri_b - 1.0)) *
           (-(2.0 * alpha_h * Ri_b - 1.0) - 
           pow(1.0 + 4.0 * (alpha_h - alpha_m) * Ri_b / Pr_0, 0.5));
    psi_m = -alpha_m * zeta * (1.0 - z_0 / z_atm);
    psi_h = -alpha_h * zeta * (1.0 - z_0 / z_atm);
  
  } else {
    
    sb = Ri_b / Pr_0;
    Qb = 1.0 / 9.0 * (1.0 / (beta_m * beta_m) + 
      3.0 * beta_h / beta_m * sb * sb);
    Pb = 1.0 / 54.0 * (-2.0 / pow(beta_m, 3.0) + 
      9 / beta_m * (-beta_h / beta_m + 3.0) * sb * sb);
    thetab = acos(Pb / pow(Qb, 1.5));
    Tb = pow(sqrt(pow(Pb, 2.0) + pow(Qb, 3.0)) + abs(Pb), 1.0/3.0);
    
    if (Ri_b > Ri_c) {
      zeta = z_atm / (z_atm - z_0) * log(z_atm / z_0) *
        (-(Tb + Qb / Tb) + 1 / (3.0 * beta_m));
    } else {
      zeta = z_atm / (z_atm - z_0) * log(z_atm / z_0) * 
        (-2.0 * sqrt(Qb) * cos(thetab / 3.0) + 1.0 / (3.0 * beta_m));
    }

    x = pow(1.0 - beta_m * zeta, 0.25);
    x_0 = pow(1.0 - beta_m * zeta * z_0 / z_atm, 0.25);
    y = pow(1.0 - beta_h * zeta, 0.5);
    y_0 = pow(1.0 - beta_h * zeta * z_0 / z_atm, 0.5);
    psi_m = 2.0 * log((1.0 + x) / (1.0 + x_0)) +
      log((1.0 + x * x) / (1.0 + x_0 * x_0)) -
      2.0 * atan(x) + 2.0 * atan(x_0);
    psi_h = 2.0 * log((1.0 + y) / (1.0 + y_0));
  }

  // Calculate surface stresses
  double tau_on_u = -k_vK * k_vK * pow(log(z_atm / z_0) - psi_m, -2.0) * u_mag;
  *taux = tau_on_u * u_atm;
  *tauy = tau_on_u * v_atm;
  
  // Calculate aerodynamic resistance
  double r_a = Pr_0 * (log(z_atm / z_0) - psi_m) * (log(z_atm / z_0) - psi_h) / 
    (k_vK * k_vK * u_mag);
  
  // Calculate soil resistance
  double r_s = r_sfc * (phi_fc - phi_pwp) / (phi - phi_pwp);

  // Calculate sensible heat flux
  *shf = rho * c_p / r_a * (theta_s - theta_atm);

  // Calculate latent heat flux
  *lhf = rho * L_v / (r_a + r_s) * (qstar_s - q_atm);

}

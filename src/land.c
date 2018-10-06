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
#define ABS(X) ((X) > 0 ? (X) : (-(X)))

/* Private constants */
/** Neutral stability Prandtl number */
static const REAL Pr_0 = 0.74;
/** von Karman constant */
static const REAL k_vK = 0.40;
/** Profile function parameters */
static const REAL alpha_m = 4.7;
static const REAL beta_m = 15.0;
static const REAL alpha_h = 6.35;
static const REAL beta_h = 9.0;
/** Maximum bulk Richardson number */
static const REAL Ri_b_max = 0.21;
/** Minimum surface-atmosphere disequilibria */
static const REAL dtheta_min = 1e-6;
static const REAL dq_min = 1e-6;

/**
 * Implements land_flux
 */
void simple_land_flux(REAL grav,
                      REAL theta_s, REAL theta_atm,
                      REAL qstar_s, REAL q_atm,
                      REAL u_atm, REAL v_atm,
                      REAL u_min,
                      REAL phi, REAL phi_fc,
                      REAL phi_pwp, REAL r_sfc,
	              REAL z_atm, REAL z_0,
	              REAL *Ri_b_out, REAL *zeta_out, 
	              REAL *C_k_out, REAL *C_q_out,
                      REAL *C_d_out,
                      REAL *shf_out, REAL *lhf_out,
                      REAL *taux_out, REAL *tauy_out) {

  // Check disequilibria for near-zero values
  REAL dtheta = theta_atm - theta_s;
  if (ABS(dtheta) < dtheta_min) {
    dtheta = dtheta_min;
  }
  REAL dq = q_atm - qstar_s;
  if (ABS(dq) < dq_min) {
    dq = dq_min;
  }

  // Calculate bulk Richardson number
  // Includes a check for Ri_b > alpha_m
  REAL u_mag = sqrt(u_atm*u_atm + v_atm*v_atm);
  if (u_mag < u_min) { 
    u_mag = u_min; 
  }
  REAL Ri_b = grav / theta_s * dtheta * z_atm / (u_mag * u_mag);
  if (Ri_b > Ri_b_max) {
    Ri_b = Ri_b_max;
  }

  // Calculate M-O stability parameter
  REAL zeta;
  REAL sb;
  REAL Qb;
  REAL Pb;
  REAL thetab;
  REAL Tb;
  REAL psi_m;
  REAL psi_h;
  REAL x;
  REAL x_0;
  REAL y;
  REAL y_0;
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
    Tb = pow(sqrt(pow(Pb, 2.0) - pow(Qb, 3.0)) + ABS(Pb), 1.0/3.0);
    
    if (pow(Qb, 3.0) - pow(Pb, 2.0) < 0) {
      zeta = z_atm / (z_atm - z_0) * log(z_atm / z_0) *
        (-(Tb + Qb / Tb) + 1.0 / (3.0 * beta_m));
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

  // Calculate exchange coefficients
  REAL u_star = k_vK * u_mag / (log(z_atm / z_0) - psi_m);
  REAL theta_star = k_vK * dtheta / Pr_0 / (log(z_atm / z_0) - psi_h);
  REAL C_d = u_star * u_star / (u_mag * u_mag);
  REAL C_k = u_star * theta_star / (u_mag * dtheta);

  // Calculate surface stresses
  REAL taux = -C_d * u_mag * u_atm;
  REAL tauy = -C_d * u_mag * v_atm;
  
  // Calculate water vapor exchange coefficient including soil resistance
  REAL C_q =
    C_k * (phi - phi_pwp) / 
    (phi - phi_pwp + C_k * u_mag * r_sfc * (phi_fc - phi_pwp));

  // Calculate sensible heat flux
  REAL shf = -C_k * u_mag * dtheta;

  // Calculate latent heat flux
  REAL lhf = -C_q * u_mag * dq;

  // Set outputs
  *Ri_b_out = Ri_b;
  *zeta_out = zeta;
  *C_d_out = C_d;
  *C_k_out = C_k;
  *C_q_out = C_q;
  *shf_out = shf;
  *lhf_out = lhf;
  *taux_out = taux;
  *tauy_out = tauy;

}

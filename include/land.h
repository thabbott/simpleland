#ifndef _LAND_H

#define land_flux simple_land_flux

/**
 * Floating point type definition
 */
typedef double REAL;


/**
 * land_flux
 *
 * Compute surface momentum, sensible heat, and latent heat fluxes. All
 * parameters should be given in MKS units.
 *
 * @param grav gravitational acceleration
 * @param theta_s potential temperature of surface
 * @param theta_atm potential temperature at lowest model level
 * @param qstar_s saturation specific humidity at temperature of surface
 * @param q_atm specific humidity at lowest model level
 * @param u_atm u velocity at lowest model level
 * @param v_atm v velocity at lowest model level
 * @param u_min minimum wind speed for surface flux computation
 * @param phi surface soil moisture
 * @param phi_fc surface soil moisture at field capacity
 * @param ph_pwp surface soil moisture at permanent wilting point
 * @param r_sfc soil resistance at field capacity
 * @param z_atm distance between surface and lowest model level
 * @param z_0 surface roughness length
 * @param zeta_out Monin-Obukhov stability parameter (output, nondim)
 * @param C_k_out Surface exchange coefficient (output, nondim)
 * @param C_q_out Surface exchange coefficient for water vapor (output, nondim)
 * @param C_d_out Surface drag coefficient (output, nondim)
 * @param shf_out surface sensible heat flux (output, K m/s)
 * @param lhf_out surface latent heat flux (output, kg/kg m/s)
 * @param taux_out surface stress in u direction (output, m^2/s^2)
 * @param tauy_out surface stress in v direction (output, m^2/s^2)
 */
extern void land_flux(REAL grav,
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
                      REAL *taux_out, REAL *tauy_out);

#define _LAND_H
#endif

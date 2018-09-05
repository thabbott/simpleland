#ifndef _LAND_H

#define land_flux simple_land_flux

/**
 * land_flux
 *
 * Compute surface momentum, sensible heat, and latent heat fluxes. All
 * parameters should be given in MKS units.
 *
 * @param rho air density at lowest model level
 * @param grav gravitational acceleration
 * @param c_p atmospheric heat capacity at constant pressure
 * @param L_v latent heat of vaporization
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
 * @param C_d_out Surface drag coefficient (output, nondim)
 * @param shf_out surface sensible heat flux (output, W/m^2)
 * @param lhf_out surface latent heat flux (output, W/m^2)
 * @param taux_out surface stress in u direction (output, m^2/s^2)
 * @param tauy_out surface stress in v direction (output, m^2/s^2)
 */
extern void land_flux(double rho, double grav,
                      double c_p, double L_v,
                      double theta_s, double theta_atm,
                      double qstar_s, double q_atm,
                      double u_atm, double v_atm,
                      double u_min,
                      double phi, double phi_fc,
                      double phi_pwp, double r_sfc,
		      double z_atm, double z_0,
		      double *Ri_b_out, double *zeta_out, 
		      double *C_k_out, double *C_d_out,
                      double *shf_out, double *lhf_out,
                      double *taux_out, double *tauy_out);

#define _LAND_H
#endif

cimport c_land

def land_flux(args):
    """
    Compute surface fluxes of momentum, sensible heat, and latent heat

    Inputs
    ------
    args: dict
        must contain keys for all input arguments listed in land.h:land_flux

    Outputs
    -------
    shf: scalar
        surface sensible heat flux
    lhf: scalar
        surface latent heat flx
    taux: scalar
        surface stress in u direction
    tauy: scalar
        surface stress in v direction
    """

    cdef double shf
    cdef double lhf
    cdef double taux
    cdef double tauy

    cdef double rho = args['rho']
    cdef double grav = args['grav']
    cdef double c_p = args['c_p']
    cdef double L_v = args['L_v']
    cdef double theta_s = args['theta_s']
    cdef double theta_atm = args['theta_atm']
    cdef double qstar_s = args['qstar_s']
    cdef double q_atm = args['q_atm']
    cdef double u_atm = args['u_atm']
    cdef double v_atm = args['v_atm']
    cdef double u_min = args['u_min']
    cdef double phi = args['u_phi']
    cdef double phi_fc = args['phi_fc']
    cdef double phi_pwp = args['phi_pwp']
    cdef double r_sfc = args['r_sfc']
    cdef double z_atm = args['z_atm']
    cdef double z_0 = args['z_0']

    c_land.land_flux(rho, grav, c_p, L_v, theta_s, theta_atm,
              qstar_s, q_atm, u_atm, v_atm, u_min,
              phi, phi_fc, phi_pwp, r_sfc, z_atm, z_0,
              &shf, &lhf, &taux, &tauy)

    return shf, lhf, taux, tauy

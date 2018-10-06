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
        surface sensible heat flux (W/m^2)
    lhf: scalar
        surface latent heat flx (W/m^2)
    taux: scalar
        surface stress in u direction (m^2/s^2)
    tauy: scalar
        surface stress in v direction (m^2/s^2)
    C_k: scalar
        surface exchange coefficient (nondim)
    C_d: scalar
        surface drag coefficient (nondim)
    zeta: scalar
        Monin-Obukhov stability parameter (nondim)
    Ri_b: scalar
        Bulk Richardson number (nondim)
    """

    cdef double Ri_b
    cdef double zeta
    cdef double C_k
    cdef double C_q
    cdef double C_d
    cdef double shf
    cdef double lhf
    cdef double taux
    cdef double tauy

    cdef double grav = args['grav']
    cdef double theta_s = args['theta_s']
    cdef double theta_atm = args['theta_atm']
    cdef double qstar_s = args['qstar_s']
    cdef double q_atm = args['q_atm']
    cdef double u_atm = args['u_atm']
    cdef double v_atm = args['v_atm']
    cdef double u_min = args['u_min']
    cdef double phi = args['phi']
    cdef double phi_fc = args['phi_fc']
    cdef double phi_pwp = args['phi_pwp']
    cdef double r_sfc = args['r_sfc']
    cdef double z_atm = args['z_atm']
    cdef double z_0 = args['z_0']
    cdef double zeta_max = args['zeta_max']

    c_land.land_flux(grav, theta_s, theta_atm,
              qstar_s, q_atm, u_atm, v_atm, u_min,
              phi, phi_fc, phi_pwp, r_sfc, z_atm, z_0, zeta_max,
              &Ri_b, &zeta, &C_k, &C_q, &C_d, &shf, &lhf, &taux, &tauy)

    return shf, lhf, taux, tauy, C_k, C_d, zeta, Ri_b

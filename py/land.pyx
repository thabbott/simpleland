cdef extern from "land.h":
    void land_flux(double rho, double grav,
                   double c_p, double L_v,
                   double theta_s, double theta_atm,
                   double qstar_s, double q_atm,
                   double u_atm, double v_atm,
                   double u_min,
                   double phi, double phi_fc,
                   double phi_pwp, double r_sfc,
                   double z_atm, double z_0,
                   double *shf, double *lhf,
                   double *taux, double *tauy);

cpdef py_land_flux(rho, grav,
                 c_p, L_v,
                 theta_s, theta_atm,
                 qstar_s, q_atm,
                 u_atm, v_atm,
                 u_min,
                 phi, phi_fc,
                 phi_pwp, r_sfc,
                 z_atm, z_0):
    cdef double shf
    cdef double lhf
    cdef double taux
    cdef double tauy

    land_flux(rho, grav, c_p, L_v, theta_s, theta_atm,
              qstar_s, q_atm, u_atm, v_atm, u_min,
              phi, phi_fc, phi_pwp, r_sfc, z_atm, z_0,
              &shf, &lhf, &taux, &tauy)

    return shf, lhf, taux, tauy

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

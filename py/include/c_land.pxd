cdef extern from "land.h":
    void land_flux(double grav,
                   double theta_s, double theta_atm,
                   double qstar_s, double q_atm,
                   double u_atm, double v_atm,
                   double u_min,
                   double phi, double phi_fc,
                   double phi_pwp, double r_sfc,
                   double z_atm, double z_0, double zeta_max,
                   double *Ri_b_out, double *zeta_out,
                   double *C_k_out, double *C_q_out,
                   double *C_d_out,
                   double *shf_out, double *lhf_out,
                   double *taux_out, double *tauy_out)

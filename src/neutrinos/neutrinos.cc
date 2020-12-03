#include "gadgetconfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

#include "pcu.h"
#include "neutrinos.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../pm/pm_periodic.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

double nulinear::tau_t_eV(int t) {
    if (N_tau==0) return 0.0;
    static int init = 0;
    static double *tau_table_eV;
    
    if(!init) {
        tau_table_eV = (double *)Mem.mymalloc_movable(&tau_table_eV, "tau_table_eV", N_tau * sizeof(double));
        gsl_interp_accel *spline_accel = gsl_interp_accel_alloc();
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen,pcu_N);
        gsl_spline_init(spline,pcu_prob,pcu_tau,pcu_N);
        
        for(int t=0; t<N_tau; t++){
            double prob = (0.5+t) / N_tau;
            tau_table_eV[t] = gsl_spline_eval(spline,prob,spline_accel);
        }
        
        gsl_spline_free(spline);
        gsl_interp_accel_free(spline_accel);
        init = 1;
    }
    
    if(t == -1){
        Mem.myfree(tau_table_eV);
        init = 0;
        return 0;
    }
    return tau_table_eV[t];
}

double nulinear::v2_t_eta_REL(int t, double eta) {
    double m_aeta_tau = m_nu_eV * aeta_in*exp(eta) / Nulinear.tau_t_eV(t);
    return 1.0 / (1.0 + m_aeta_tau*m_aeta_tau);
}

double nulinear::Ec_t_eta(int t, double eta) { return 1.0 / (aeta_in*exp(eta) ); }

double nulinear::Ec_t_eta_REL(int t, double eta) {
    double vt2 = Nulinear.v2_t_eta_REL(t,eta), aeta = aeta_in*exp(eta);
    if(1-vt2 < 1e-12){
        double ma_tau = m_nu_eV * aeta / Nulinear.tau_t_eV(t);
        return sqrt(1.0 + ma_tau*ma_tau) / (aeta*ma_tau);
    }
    return 1.0 / (aeta * sqrt(1.0 - vt2));
}

double nulinear::eta_convert(double a) {
    return log(a/aeta_in);
}

double nulinear::Ec_de_eta(double eta) {
    double aeta = aeta_in * exp(eta);
    return pow(aeta, -1.0 -3.0*(w0_eos_de + wa_eos_de)) * exp(3.0*wa_eos_de*(aeta-1.0));
}

double nulinear::Hc2_Hc02_eta(double eta) {
    //scale factor
    double aeta = aeta_in * exp(eta), aeta2 = aeta*aeta, Ec_de = Nulinear.Ec_de_eta(eta);
    
    //sum Omega_{t,0} aeta^2 rho_t(eta)/rho_t_0 over CDM, photons, and DE
    double sum_OEc = All.Omega0/aeta + Omega_rel_0/aeta2 + All.OmegaLambda*Ec_de;
    
    //neutrinos
    for(int t=0; t<N_tau; t++) sum_OEc += Omega_nu_t_0 * Nulinear.Ec_t_eta_REL(t,eta);
    
    return sum_OEc;
}

double nulinear::Hubble_Extern (double a) {
    return sqrt(Nulinear.Hc2_Hc02_eta(Nulinear.eta_convert(a)))/a;
}

double nulinear::OF_eta(int F, double eta) {
    double Hc02_Hc2 = 1.0 / Nulinear.Hc2_Hc02_eta(eta), aeta = aeta_in*exp(eta);

    if(F == N_tau) // CDM
      return All.Omega0 * pow(aeta, -1.0-3.0*w_eos_cdm) * Hc02_Hc2;
    else if(F == N_tau+1) // photon + massless nu
      return Omega_rel_0 * pow(aeta, -1.0-3.0*w_eos_gam) * Hc02_Hc2;
    else if(F == N_tau+2) // dark energy, assumed Lambda
      return Omega_de_0 * Nulinear.Ec_de_eta(eta) * Hc02_Hc2;
    else if(F<0 || F>N_tau+2) return 0.0; // no fluids  should have these indices
    return Omega_nu_t_0 * Nulinear.Ec_t_eta(F,eta) * Hc02_Hc2;
}

double nulinear::Hc_eta(double eta) { return Hc0h * sqrt(Nulinear.Hc2_Hc02_eta(eta)); }

double nulinear::poisson_mod_fac(double k, double a, double dcb) {
    double fnu = All.OmegaNu/(All.Omega0+All.OmegaNu);
    double kfs = 1.5 * sqrt(All.Time * (All.Omega0 + All.OmegaNu)) * Nulinear.m_nu_eV_parser();

    return ((k + kfs) * (k + kfs)) / ((k + kfs) * (k + kfs) - kfs * kfs * fnu); 
}

double nulinear::m_nu_eV_parser(void) {
    return m_nu_eV;
}

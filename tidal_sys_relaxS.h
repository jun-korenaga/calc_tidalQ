/*
 * tidal_sys_relaxS.h (banded radial solver)
 *
 * note: this is specific to the Earth-Moon system
 *       with a focus on tidal dissipation within Earth.
 *
 * Jun Korenaga
 * January 2025
 */

#ifndef _TIDAL_SYS_S_H_
#define _TIDAL_SYS_S_H_

#include "const_tidal.h"
#include "earth_propS.h"
#include <iostream>

class TidalSystem {
public:
    // for a circular orbit
    TidalSystem(EarthProperties& _earth);
    int calcLoveNumbers(double omega_t, int order,
			complex<double>& h2, complex<double>& k2,
			double& delta_h, double& delta_k,
			double& tidalQ_h, double& tidalQ_k,
			double& error_val);
    void dumpSolutions(const char* fn_root);
    void inertiaOff(){ inertia_on = false; inertia_CMB = false; }
    void noInertiaCMB(){ inertia_CMB = false; }

private:
    double cabs(complex<double>);
    Array1d< complex<double> > set_initial_inertia(int k, double omegaa_t);
    double sbf(int n, double x);
    double isbf_ex(int n, double x);
    double doublefactorial(int n);
    
    double solveRadialSystem();

    bool inertia_on, inertia_CMB, dump_solution;
    const char* dump_outroot;
    int nr;
    int dump_nfile;
    complex<double> cmplx_zero;
    double omega_ep, omega_mp, dist_p, I_e, L;
    double rhoc, g_r, rho_r, Phi_r, mu_r, degL;
    double max_tidalQ;
    double tidal_BC_n;
    complex<double> CMB_a0, CMB_b0, CMB_c0, CMB_d0;
    complex<double> CMB_a, CMB_b, CMB_c, CMB_d;

    EarthProperties& earth;
    Array1d<double> r, r_n, rho0, rho0_n, g_n, kappa, kappa_n;
    Array1d< complex<double> > lambda, lambda_n, mu, mu_n, beta_n;
    Array1d< Array1d< complex<double> > > yref;
    Array2d< complex<double> > matI;
    Array2d< Array1d< complex<double> > > A;
};

#endif /* _TIDAL_SYS_S_H_ */

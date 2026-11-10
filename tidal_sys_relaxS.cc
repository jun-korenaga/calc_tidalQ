/*
 * tidal_sys_relaxS.cc
 *
 * Jun Korenaga
 * January 2025
 */

#include "tidal_sys_relaxS.h"
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

// Use LAPACKE's C++ complex type overrides for its array arguments.
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>

TidalSystem::TidalSystem(EarthProperties& _earth)
    : earth(_earth)
{
    // default values
    inertia_on = true;
    inertia_CMB = true;
    dump_solution = false;
    max_tidalQ = 1e10;
    
    // present-day state
    omega_ep = 2*tidal_pi/day2sec; // Earth angular rotation rate [rad/s]
    omega_mp = 2*tidal_pi/(27.322*day2sec); // Moon angular rotation rate
    dist_p = pow(G*(M_e+M_m)/pow(omega_mp,2.0),1.0/3.0); // Earth-Moon distance

    // time invariant parameters
    I_e = 1.0/3.0*M_e*RE*RE; // Earth's moment of inertia
    L = I_e*omega_ep+M_m*sqrt(G*(M_e+M_m)*dist_p); // total angular momentum

    earth.setRadius(r);
    earth.setDensity(rho0, rhoc);
    earth.setGravity(g_n);
    earth.setElasticModuli(kappa, lambda, mu);
    nr = r.size();
	
    g_r = g_n.back(); // reference gravity
    rho_r = M_e/(4.0*tidal_pi/3.0*pow(RE,3.0)); // reference density    
    Phi_r = g_r*RE; // potential scale
    mu_r = rho_r*g_r*RE; // stress scale

    r_n = r; r_n /= RE;
    rho0_n = rho0; rho0_n /= rho_r;
    g_n /= g_r;
    kappa_n = kappa; kappa_n /= mu_r;

    degL = 2.0; // degree of spherical harmonics

    cmplx_zero.real(0.0);
    cmplx_zero.imag(0.0);
    A.resize(6,6);
    for (int i=1; i<=6; i++){
	for (int j=1; j<=6; j++){
	    A(i,j).resize(nr); A(i,j) = cmplx_zero;
	}
    }

    // static CMB boundary condition prep. (normalized)
    matI.resize(6,3);
    matI = cmplx_zero;
    double rc = r_n(1)*RE;
    double gc = g_n(1)*g_r;
    double rcL = pow(rc,degL);
    matI(1,1).real(-rcL/(gc*RE));
    matI(1,3).real(1.0/RE);
    matI(2,2) = matI(1,3);
    matI(3,3).real(gc*rhoc/mu_r);
    matI(5,1).real(rcL/Phi_r);
    matI(6,1).real(((2.0*degL+1.0)*pow(rc,(degL-1.0))
		    -4.0*tidal_pi*G*rhoc/gc*rcL)/g_r);
    matI(6,3).real(4.0*tidal_pi*G*rhoc/g_r);

    yref.resize(nr);
    for (int i=1; i<=nr; i++){
	yref(i).resize(6);
    }
}

void TidalSystem::dumpSolutions(const char* fn_root)
{
    dump_solution = true;
    dump_outroot = fn_root;
    dump_nfile = 0;
}

int TidalSystem::calcLoveNumbers(double omega_t, int order,
				 complex<double>& h2,
				 complex<double>& k2,
				 double& delta_h, double& delta_k,
				 double& tidalQ_h, double& tidalQ_k,
				 double& err_val)
{
    lambda_n = lambda; lambda_n /= mu_r;
    mu_n = mu; mu_n /= mu_r;
    beta_n = mu_n; beta_n *= 2.0; beta_n += lambda_n;
    double omega_t2 = omega_t*omega_t;
    double RE2 = RE*RE;

    for (int i=1; i<=nr; i++){
	double inertia_term = 0.0;
	if (inertia_on) inertia_term =  (RE2/mu_r)*rho0(i)*omega_t2;

	A(1,1)(i) = -2.0*lambda_n(i)/(r_n(i)*beta_n(i));
	A(1,2)(i) =  degL*(degL+1.0)*lambda_n(i)/(r_n(i)*beta_n(i));
	A(1,3)(i) =  1.0/beta_n(i);

	A(2,1)(i) =  -1.0/r_n(i);
	A(2,2)(i) =  1.0/r_n(i);
	A(2,4)(i) =  1.0/mu_n(i);

	A(3,1)(i) =  4.0/r_n(i)*(3.0*kappa_n(i)*mu_n(i)/(r_n(i)*beta_n(i))
				 -rho0_n(i)*g_n(i)) - inertia_term;
	A(3,2)(i) =  degL*(degL+1.0)/r_n(i)*(rho0_n(i)*g_n(i)
					     -6.0*kappa_n(i)*mu_n(i)/(r_n(i)*beta_n(i)));
	A(3,3)(i) =  -4.0*mu_n(i)/(r_n(i)*beta_n(i));
	A(3,4)(i) =  degL*(degL+1.0)/r_n(i);
	A(3,5)(i) =  -rho0_n(i)*(degL+1.0)/r_n(i);
	A(3,6)(i) =  rho0_n(i);

	A(4,1)(i) =  1.0/r_n(i)*(rho0_n(i)*g_n(i)
				 -6.0*kappa_n(i)*mu_n(i)/(r_n(i)*beta_n(i)));
	A(4,2)(i) =  2.0*mu_n(i)/(r_n(i)*r_n(i))*(degL*(degL+1.0)
						  *(1.0+lambda_n(i)/beta_n(i))-1.0)
	    - inertia_term;
	A(4,3)(i) =  -lambda_n(i)/(r_n(i)*beta_n(i));
	A(4,4)(i) =  -3.0/r_n(i);
	A(4,5)(i) =  rho0_n(i)/r_n(i);
	    
	A(5,1)(i) =  -3.0*rho0_n(i);
	A(5,5)(i) =  -(degL+1.0)/r_n(i);
	A(5,6)(i) = 1.0;

	A(6,1)(i) =  -3.0*rho0_n(i)*(degL+1.0)/r_n(i);
	A(6,2)(i) =  A(6,1)(i)*(-degL);
	A(6,6)(i) =  (degL-1.0)/r_n(i);
    }

    // for CMB boundary conditions
    if (inertia_on && inertia_CMB){
	Array1d< Array1d< complex<double> > > ivec(2);
	for (int i=1; i<=2; i++){
	    ivec(i).resize(6);
	    ivec(i) = set_initial_inertia(i,omega_t);
	    double norm=0.0;
	    for (int I=1; I<=6; I++){
		double scale = cabs(ivec(i)(I));
		norm += scale*scale;
	    }
	    norm = sqrt(norm);
	    ivec(i) /= norm;
	}
	CMB_a = ivec(1)(5)/ivec(1)(1);
	CMB_b = ivec(2)(5)/ivec(2)(3)
	    - ivec(2)(1)*ivec(1)(5)/(ivec(1)(1)*ivec(2)(3));
	CMB_c = ivec(1)(6)/ivec(1)(1);
	CMB_d = ivec(2)(6)/ivec(2)(3)
	    - ivec(2)(1)*ivec(1)(6)/(ivec(1)(1)*ivec(2)(3));
    }else{ // static BC
	CMB_a0 = matI(1,3)/matI(3,3);
	CMB_b0 = matI(1,1)/matI(5,1);
	CMB_c0 = matI(6,3)/matI(3,3);
	CMB_d0 = matI(6,1)/matI(5,1);
    }

    double dist = dist_p; // note: this doesn't matter
    double gamma = M_m/M_e*pow(RE/dist,3.0)*RE;
    double pfac = 0;
    if (order==2){
	pfac = 1.0/4.0;
    }else if (order==1){
	pfac = 1.0/2.0;
    }else if (order==0){
	pfac = 1.0;
    }else{
	cerr << "unimplemented order (" << order << ") detected.\n";
	exit(1);
    }

    double tidal_potential = -pfac*gamma*g_r;
    double tidal_BC = (2.0*degL+1.0)/RE*tidal_potential;
    tidal_BC_n = tidal_BC/g_r;

    err_val = solveRadialSystem();

    if (dump_solution){
	dump_nfile++;
        std::ostringstream filename;
        filename << dump_outroot << '.' << std::setw(2) << std::setfill('0')
                 << dump_nfile << ".dat";
        const std::string path = filename.str();
        std::ofstream dump(path);
        if (!dump) throw std::runtime_error("Cannot open solution dump: " + path);

        for (int i=1; i<=nr; i++){
            dump << r(i) << " "
                 << yref(i)(1).real()*RE << " " << yref(i)(1).imag()*RE << " "
                 << yref(i)(2).real()*RE << " " << yref(i)(2).imag()*RE << " "
                 << yref(i)(3).real()*mu_r << " " << yref(i)(3).imag()*mu_r << " "
                 << yref(i)(4).real()*mu_r << " " << yref(i)(4).imag()*mu_r << " "
                 << yref(i)(5).real()*Phi_r << " " << yref(i)(5).imag()*Phi_r << " "
                 << yref(i)(6).real()*g_r << " " << yref(i)(6).imag()*g_r << " "
                 << mu(i).real() << " " << mu(i).imag() << '\n';
        }
        dump.close();
        if (!dump) throw std::runtime_error("Cannot write solution dump: " + path);
    }
    
    //
    // calculate Love numbers and tidal Q
    //
    h2 = yref(nr)(1)*RE/(pfac*gamma);
    k2 = (-yref(nr)(5)*Phi_r)/(pfac*gamma*g_r)-1.0;
    delta_h = atan2(-h2.imag(),h2.real());
    delta_k = atan2(-k2.imag(),k2.real());

    tidalQ_h = max_tidalQ;
    tidalQ_k = max_tidalQ;
    if (abs(h2.imag())>1e-10) tidalQ_h = abs(h2.real())/abs(h2.imag());
    if (abs(k2.imag())>1e-10) tidalQ_k = abs(k2.real())/abs(k2.imag());

    return 1; // One direct solve.
}

double TidalSystem::solveRadialSystem()
{
    constexpr int components = 6;
    constexpr int core_conditions = 3;
    // Each interior row couples two neighboring nodes. With three core
    // rows first, the farthest column is eight places from the diagonal.
    constexpr int bandwidth = components + core_conditions - 1;
    constexpr int leading_dim = 3 * bandwidth + 1;
    const int size = components * nr;
    using Complex = std::complex<double>;
    std::vector<Complex> band(static_cast<std::size_t>(size) * leading_dim);
    std::vector<Complex> forcing(size);

    // LAPACK stores diagonals in rows, leaving room above for pivot fill-in.
    const auto coefficient = [&](int row, int column) -> Complex& {
        return band[static_cast<std::size_t>(column) * leading_dim +
                    (2 * bandwidth + (row - column))];
    };

    // Core–mantle boundary at the first node.
    if (inertia_on && inertia_CMB) {
        coefficient(0, 0) = -CMB_a;
        coefficient(0, 2) = -CMB_b;
        coefficient(0, 4) = 1.0;
        coefficient(2, 0) = -CMB_c;
        coefficient(2, 2) = -CMB_d;
        coefficient(2, 5) = 1.0;
    } else {
        coefficient(0, 0) = 1.0;
        coefficient(0, 2) = -CMB_a0;
        coefficient(0, 4) = -CMB_b0;
        coefficient(2, 2) = -CMB_c0;
        coefficient(2, 4) = -CMB_d0;
        coefficient(2, 5) = 1.0;
    }
    coefficient(1, 3) = 1.0;

    // Interior equations use the original centered discretization:
    // y_right - y_left = (dr / 2) * mean(A_left, A_right) * (y_left + y_right).
    for (int node = 2; node <= nr; ++node) {
        const int left = components * (node - 2);
        const int right = left + components;
        const double half_spacing = 0.5 * (r_n(node) - r_n(node - 1));
        for (int eq = 0; eq < components; ++eq) {
            const int row = left + core_conditions + eq;
            for (int var = 0; var < components; ++var) {
                const Complex average =
                    0.5 * (A(eq + 1, var + 1)(node - 1) + A(eq + 1, var + 1)(node));
                const double identity = eq == var ? 1.0 : 0.0;
                coefficient(row, left + var) = -identity - half_spacing * average;
                coefficient(row, right + var) = identity - half_spacing * average;
            }
        }
    }

    // Free surface tractions and the imposed tide-raising potential.
    const int surface = size - components;
    coefficient(surface + 3, surface + 2) = 1.0;
    coefficient(surface + 4, surface + 3) = 1.0;
    coefficient(surface + 5, surface + 5) = 1.0;
    forcing[surface + 5] = tidal_BC_n;

    const auto finite = [](Complex value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
    };
    if (!std::all_of(band.begin(), band.end(), finite) ||
        !finite(forcing.back())) {
        throw std::runtime_error("Nonfinite coefficient in the radial tidal system");
    }

    // ZGBSV overwrites both inputs. Keep the original system for the residual.
    auto factors = band;
    auto solution = forcing;
    std::vector<lapack_int> pivots(size);
    const lapack_int info = LAPACKE_zgbsv(
        LAPACK_COL_MAJOR, size, bandwidth, bandwidth, 1, factors.data(),
        leading_dim, pivots.data(), solution.data(), size);
    if (info > 0) {
        throw std::runtime_error("Singular radial tidal system at pivot " +
                                 std::to_string(info));
    }
    if (info < 0) {
        throw std::runtime_error("LAPACKE_zgbsv failed with status " +
                                 std::to_string(info));
    }
    if (!std::all_of(solution.begin(), solution.end(), finite)) {
        throw std::runtime_error("Nonfinite radial tidal solution");
    }

    double residual_norm = 0.0;
    double matrix_norm = 0.0;
    double solution_norm = 0.0;
    for (int row = 0; row < size; ++row) {
        Complex residual = -forcing[row];
        double row_norm = 0.0;
        const int last = row + std::min(bandwidth, size - 1 - row);
        for (int column = std::max(0, row - bandwidth); column <= last; ++column) {
            const Complex value = coefficient(row, column);
            residual += value * solution[column];
            row_norm += std::abs(value);
        }
        const double magnitude = std::abs(residual);
        if (!std::isfinite(magnitude)) {
            throw std::runtime_error("Nonfinite radial tidal residual");
        }
        residual_norm = std::max(residual_norm, magnitude);
        matrix_norm = std::max(matrix_norm, row_norm);
        solution_norm = std::max(solution_norm, std::abs(solution[row]));
        yref(row / components + 1)(row % components + 1) = solution[row];
    }

    const double scale = matrix_norm * solution_norm + std::abs(forcing.back());
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::runtime_error("Invalid scale for the radial tidal residual");
    }
    return residual_norm / scale;
}

Array1d< complex<double> > TidalSystem::set_initial_inertia(int kk,
							    double omega_t)
{
    Array1d< complex<double> > yCMB(6);
    yCMB = cmplx_zero;

    double rc = r_n(1)*RE;
    double alpha = 8060.00; // PREM P-wave velocity just below CMB

    double omega2 = omega_t*omega_t;
    double gammaCMB = 4.0*tidal_pi*G*rhoc/3.0;

    if (kk==1){
	double r_l_1 = pow(rc,double(degL-1));
	yCMB(1) = degL*r_l_1;
	yCMB(2) = r_l_1;
	yCMB(5) = -(degL*gammaCMB-omega2)*r_l_1*rc;
	// cerr << "y5 calc " << degL << " " << gammaCMB << " " << omega2 << endl;
	yCMB(6) = -(2.0*(degL-1.0)*degL*gammaCMB-(2.0*degL+1.0)*omega2)*r_l_1;
    }else if (kk==2){
	double alpha2 = alpha*alpha;
	complex<double> ksq = 1.0/alpha2*(omega2+4.0*gammaCMB
					  -degL*(degL+1)*gammaCMB*gammaCMB/omega2);
	complex<double> k = sqrt(ksq);
	double f = -omega2/gammaCMB;
	double h = f-(degL+1.0);
	double r_l = pow(rc,double(degL));
	double r_lp1 = r_l*rc;
	double r_lp2 = r_lp1*rc;
	double lambda_c = alpha2*rhoc; // alpha = sqrt(lambda/rho)
	double rc2 = rc*rc;

	if (k.imag()==0){
	    double xc = k.real()*rc;
	    double phi_l = doublefactorial(2*degL+1)
		*pow(xc,-double(degL))*sbf(degL,xc);
	    double phi_l1 = doublefactorial(2*(degL+1)+1)
		*pow(xc,-double(degL+1))*sbf(degL+1,xc);
	    double psi_l = 2.0*(2.0*degL+3.0)/(xc*xc)*(1-phi_l);
	    
	    yCMB(1) = -r_lp1/(2.0*degL+3.0)*(0.5*degL*h*psi_l+f*phi_l1);
	    yCMB(2) = -r_lp1/(2.0*degL+3.0)*(0.5*h*psi_l-phi_l1);
	    yCMB(3) = -lambda_c*r_l*f*phi_l;
	    yCMB(5) = -r_lp2*(alpha2*f/rc2
			      - 3.0*gammaCMB*f/(2.0*(2.0*degL+3.0))*psi_l);
	    yCMB(6) = -r_lp1*((2.0*degL+1.0)*alpha2*f/rc2
			      -3.0*gammaCMB*((2.0*degL+1.0)*f-degL*h)*psi_l
			      /(2.0*(2.0*degL+3.0)));
	}else{
	    double ixc = k.imag()*rc;
	    // note: the signs in the following assume degL==2
	    double phi_l_ex = -doublefactorial(2*degL+1)
		*pow(ixc,-double(degL))*isbf_ex(degL,ixc);
	    double phi_l1_ex = -doublefactorial(2*(degL+1)+1)
		*pow(ixc,-double(degL+1))*isbf_ex(degL+1,ixc);
	    double psi_l_ex = -2.0*(2.0*degL+3.0)/(ixc*ixc)*(exp(-ixc)-phi_l_ex);

	    yCMB(1) = -r_lp1/(2.0*degL+3.0)*(0.5*degL*h*psi_l_ex+f*phi_l1_ex);
	    yCMB(2) = -r_lp1/(2.0*degL+3.0)*(0.5*h*psi_l_ex-phi_l1_ex);
	    yCMB(3) = -lambda_c*r_l*f*phi_l_ex;
	    yCMB(5) = r_lp2*(alpha2*f/rc2*exp(-ixc)
			     - 3.0*gammaCMB*f/(2.0*(2.0*degL+3.0))*psi_l_ex);
	    yCMB(6) = (1.0/rc)*((2.0*degL+1.0)*yCMB(5)
				+3.0*degL*gammaCMB*h*r_lp2*psi_l_ex
				/(2.0*(2.0*degL+3.0)));
	    yCMB(5) *= -1;
	    yCMB(6) *= -1;
	    // yCMB(6) = -r_lp1*((2.0*degL+1.0)*alpha2*f/rc2*exp(-ixc)
	    // 		      -3.0*gammaCMB*((2.0*degL+1.0)*f-degL*h)*psi_l_ex
	    // 		      /(2.0*(2.0*degL+3.0)));
	}
    }else if (kk==3){
	yCMB(2) = RE; // this becomes unity after nondimensionalization
    }else{
	cerr << "TidalSystem::set_initial_inertia - impossible kk detected.\n";
	exit(1);
    }

    // nondimensionalization
    yCMB(1) /= RE;
    yCMB(2) /= RE;
    yCMB(3) /= mu_r;
    yCMB(4) /= mu_r;
    yCMB(5) /= Phi_r;
    yCMB(6) /= g_r;

    return yCMB;
}

double TidalSystem::sbf(int n, double x)
{
    // spherical Bessel function of the first kind
    if (n<0){
	cerr << "TidalSystem::sbf n=" << n << " is invaid.\n";
	exit(1);
    }
    
    if (x==0){
        if (n==0){
            return 1.0;  // j_0(0) = 1
        }else{
            return 0.0;  // j_n(0) = 0 for n > 0
        }
    }

    if (abs(x) < 1e-8) {
        // For very small x, use series expansion
        double term = pow(x,double(n))/double(2*n+1);
        double sum = term;
        int k = 1;
        while (abs(term)>1e-15) {
            term *= -x*x/double(2*k*(2*n+2*k+1));
            sum += term;
            k++;
        }
        return sum;
    }

    if (n==0){
        return std::sin(x)/x;
    }else if (n==1){
        return std::sin(x)/(x*x) - std::cos(x)/x;
    }else{
        double jnm1 = std::sin(x)/(x*x) - std::cos(x)/x;  // j_1(x)
        double jnm2 = std::sin(x)/x;  // j_0(x)
        for (int m = 1; m < n; ++m){
            double jn = (2 *m+1)/x*jnm1 - jnm2;
            jnm2 = jnm1;
            jnm1 = jn;
        }
        return jnm1;
    }
}

double TidalSystem::isbf_ex(int n, double x)
{
    double val;
    double rx = 1.0/x;
    double rx2 = rx*rx;
    double rx3 = rx2*rx;
    
    if (n==2){
	// j_2(i*x)/exp(x)
	val = 0.5*((3.0*rx2 - 3.0*rx3 - 1.0*rx) 
		   + exp(-2.0*x)*(3.0*rx2 + 3.0*rx3 + 1.0*rx));
    }else if (n==3){
	// j_3(i*x)/(i*exp(x))
	double rx4 = rx3*rx;
	val = 0.5*((15.0*rx4 + 6.0*rx2 - 15.0*rx3 - 1.0*rx) 
		   - exp(-2.0*x)*(15.0*rx4 + 6.0*rx2 + 15.0*rx3 + 1.0*rx));
    }else{
	cerr << "TidalSystem::isbf_ex n=" << n << " is not implemented yet.\n";
	exit(1);
    }

    return val;
}

double TidalSystem::doublefactorial(int n)
{
    if (n==0 || n==-1) {
        return 1.0;  // By definition, 0!! = 1 and (-1)!! = 1
    }
    double result = 1.0;
    for (int i=n; i>0; i-=2) {
        result *= i;
    }
    return result;
}

double TidalSystem::cabs(complex<double> c)
{
    double vr = c.real();
    double vi = c.imag();
    return sqrt(vr*vr+vi*vi);
}

	


/*
 * earth_propS.cc
 *
 * Jun Korenaga
 * January 2025
 */

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include "util.h"
#include "earth_propS.h"
#include "const_tidal.h"

EarthProperties::EarthProperties(double dr_min, double dr_max,
                                 const char* model_fn)
{
    std::ifstream input(model_fn);
    if (!input) {
        throw std::runtime_error(std::string("Cannot open model: ") + model_fn);
    }

    Array1d<double> model_z, model_rho, model_kappa, model_muR, model_muI;
    std::string line;
    std::size_t line_number = 0;
    double previous_radius = RE;
    while (std::getline(input, line)) {
        ++line_number;
        line = line.substr(0, line.find('#'));
        std::istringstream row(line);
        row >> std::ws;
        if (row.eof()) continue;

        double radius, density, bulk, shear_real, shear_imag;
        std::string extra;
        if (!(row >> radius >> density >> bulk >> shear_real >> shear_imag) ||
            (row >> extra) || !std::isfinite(radius) || !std::isfinite(density) ||
            !std::isfinite(bulk) || !std::isfinite(shear_real) ||
            !std::isfinite(shear_imag)) {
            throw std::runtime_error(std::string(model_fn) + ":" +
                                     std::to_string(line_number) +
                                     ": expected five finite numbers");
        }
        if (radius < Rc || radius > previous_radius || density <= 0.0 ||
            bulk <= 0.0 || (shear_real == 0.0 && shear_imag == 0.0)) {
            throw std::runtime_error(std::string(model_fn) + ":" +
                                     std::to_string(line_number) +
                                     ": expected descending mantle radii, positive density "
                                     "and bulk modulus, and nonzero shear modulus");
        }
        model_z.push_back(RE - radius);
        model_rho.push_back(density);
        model_kappa.push_back(bulk);
        model_muR.push_back(shear_real);
        model_muI.push_back(shear_imag);
        previous_radius = radius;
    }
    if (!input.eof() || input.bad()) {
        throw std::runtime_error(std::string("Cannot read model: ") + model_fn);
    }
    if (model_z.size() < 2 || model_z.front() != 0.0 || previous_radius != Rc) {
        throw std::runtime_error(std::string(model_fn) +
                                 ": model must span 6371000 m to 3480000 m");
    }

    double alp = log(dr_max/dr_min);
    double tmpval = 1.0;
    if (alp>1e-10){
        tmpval = alp/(exp(alp)-1.0);
    }
    double intervals = round((RE-Rc)/dr_min*tmpval);
    // The banded system has six variables per grid point.
    if (!std::isfinite(intervals) || intervals < 1 ||
        intervals >= std::numeric_limits<int>::max() / 6) {
        throw std::runtime_error("Grid spacing is out of range");
    }
    nr = static_cast<int>(intervals) + 1;
    r.resize(nr); z.resize(nr);
    for (int i=1; i<=nr; i++){
        if (alp>1e-10){
            r(i) = RE+(Rc-RE)*(exp(alp*(nr-i)/(nr-1.0))-1.0)/(exp(alp)-1.0);
        }else{
            r(i) = RE+(Rc-RE)*(nr-i)/(nr-1.0);
        }
        int j = nr-i+1;
        z(j) = RE-r(i);
    }

    // set elastic moduli etc.
    rho.resize(nr); kappa.resize(nr);
    mu.resize(nr); lambda.resize(nr); 
    for (int i=1; i<=nr; i++){
	int j = nr-i+1;
	rho(i) = intp1d(model_z,model_rho,z(j));
	kappa(i) = intp1d(model_z,model_kappa,z(j));
	mu(i).real(intp1d(model_z,model_muR,z(j)));
	mu(i).imag(intp1d(model_z,model_muI,z(j)));
	if (!std::isfinite(mu(i).real()) || !std::isfinite(mu(i).imag()) ||
            mu(i) == std::complex<double>{}) {
            throw std::runtime_error("Interpolated shear modulus is invalid at radius " +
                                     std::to_string(r(i)) + " m");
        }
        lambda(i) = kappa(i)-2.0/3.0*mu(i);
    }

    // set uniform core density
    rhoc = M_ec/(4.0*tidal_pi/3.0*pow(Rc,3.0));

    // set gravity
    g.resize(nr);
    double val = M_ec;
    double r2_prev = r(1)*r(1);
    g(1) = G*M_ec/r2_prev;
    for (int i=2; i<=nr; i++){
	double r2 = r(i)*r(i);
	double dr = r(i)-r(i-1);
	val += 2.0*tidal_pi*(r2_prev*rho(i-1)+r2*rho(i))*dr;
	g(i) = G*val/r2;
	r2_prev = r2;
    }
}

void EarthProperties::setRadius(Array1d<double>& _r)
{
    if (_r.size() != r.size()) _r.resize(r.size());
    _r = r;
}

void EarthProperties::setDepth(Array1d<double>& _z)
{
    if (_z.size() != z.size()) _z.resize(z.size());
    _z = z;
}

void EarthProperties::setGravity(Array1d<double>& _g)
{
    if (_g.size() != g.size()) _g.resize(g.size());
    _g = g;
}

void EarthProperties::setDensity(Array1d<double>& _rho, double& _rhoc)
{
    if (_rho.size() != rho.size()) _rho.resize(rho.size());
    _rho = rho;
    _rhoc = rhoc;
}

void EarthProperties::setElasticModuli(Array1d<double>& _kappa,
				       Array1d< complex<double> >& _lambda,
				       Array1d< complex<double> >& _mu)
{
    if (_kappa.size() != kappa.size()) _kappa.resize(kappa.size());
    _kappa = kappa;
    if (_lambda.size() != lambda.size()) _lambda.resize(lambda.size());
    _lambda = lambda;
    if (_mu.size() != mu.size()) _mu.resize(mu.size());
    _mu = mu;
}


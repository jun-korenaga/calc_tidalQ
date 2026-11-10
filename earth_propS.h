/*
 * earth_propS.h
 * ("S" for a simplified version)
 *
 * Jun Korenaga
 * January 2025
 */

#ifndef _EARTH_PROP_S_H_
#define _EARTH_PROP_S_H_

#include <complex>
#include "array.h"
#include "const_tidal.h"

class EarthProperties {
public:
    EarthProperties(double dr_min, double dr_max,
		    const char* model_fn);

    int nGrid(){ return nr; }
    void setRadius(Array1d<double>& _r);
    void setDepth(Array1d<double>& _z);
    void setGravity(Array1d<double>& _g);
    void setDensity(Array1d<double>& _rho, double& _rhoc);
    void setElasticModuli(Array1d<double>& _kappa,
			  Array1d< complex<double> >& _lambda,
			  Array1d< complex<double> >& _mu);

private:
    int nr;
    double rhoc;
    Array1d<double> r, rho, kappa, g;
    Array1d< complex<double> > lambda, mu;
    Array1d<double> z;
};
    

#endif /* _EARTH_PROP_S_H_ */


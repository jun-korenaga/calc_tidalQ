/*
 * util.h
 */

#ifndef _JK_UTIL_H_
#define _JK_UTIL_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <complex>
#include "array.h"

const int MaxStr = 512;

FILE *openFile(char *fn, char *mode);
int countLines(char *fn);
double intp1d(const Array1d<double>& x, const Array1d<double>& f, double xval);
double intp2d(double a, double b, double c, double d,
	      double x, double y);
complex<double> intp2d(complex<double> a, complex<double> b,
		       complex<double> c, complex<double> d,
		       double x, double y);
double cubic_spline(Array1d<double>&xs, Array1d<double>& ys,
		    Array1d<double>& zs, double x);
void check_range(double val, double val_min, double val_max,
		 string val_name, bool& gotError);
void csvd(Array2d< complex<double> >& A, Array2d< complex<double> >& U,
	  Array1d<double>& S, Array2d< complex<double> >& V);

#endif /* _JK_UTIL_H_ */

/*
 * util.cc
 */

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <ctime>
#include <cmath>
#include "util.h"
#include "error.h"

FILE *openFile(char *fn, char *mode)
{
    FILE *tmp;
    
    if ((tmp=fopen(fn, mode)) == NULL){
	fprintf(stderr, "[openFile] can't open file %s\n", fn);
	exit(1);
    }
    return tmp;
}

int countLines(char *fn)
{
    FILE *tmpfp;
    int nline=0;
    char line[MaxStr];

    tmpfp = openFile(fn, (char *) "r");
    while(fgets(line, MaxStr, tmpfp) != NULL){
	nline++;
    }
    fclose(tmpfp);

    return nline;
}

double intp1d(const Array1d<double>& x, const Array1d<double>& f, double xval)
{
    if (x.size() != f.size())
	error("intp1d::size mismatch");


    if (x.front()>x.back()){
	cerr << "intp1d requires x to be monotonically inceasing.\n";
	cerr << "x.size() = " << x.size() << ", xval = " << xval << endl;
	cerr << x.front() << " " << x.back() << " "
	     << f.front() << " " << f.back() << endl;
	exit(1);
    }

    if (xval < x.front()) return f.front();
    
    for (int i=1; i<int(x.size()); i++){
	if (xval>=x(i) && xval<x(i+1)){
	    double dx = xval-x(i);
	    double dX = x(i+1)-x(i);
	    double dF = f(i+1)-f(i);
	    return f(i)+dx/dX*dF;
	}
    }

    return f.back();
}

double intp2d(double a, double b, double c, double d,
	      double x, double y)
{
    double ab = a*(1.0-x)+b*x;
    double cd = c*(1.0-x)+d*x;
    double val = ab*(1.0-y)+cd*y;
    return val;
}

complex<double> intp2d(complex<double> a, complex<double> b,
		       complex<double> c, complex<double> d,
		       double x, double y)
{
    complex<double> ab = a*(1.0-x)+b*x;
    complex<double> cd = c*(1.0-x)+d*x;
    complex<double> val = ab*(1.0-y)+cd*y;
    return val;
}

double cubic_spline(Array1d<double>&xs, Array1d<double>& ys,
		    Array1d<double>& zs, double x)
{
    // check the range
    if (x<xs.front()) return ys.front();
    if (x>xs.back()) return ys.back();

    for (int i=2; i<=int(xs.size()); i++){
	if (x>=xs(i-1) && x<=xs(i)){
	    double h = xs(i)-xs(i-1);
	    double val = zs(i)*pow(x-xs(i-1),3.0)/(6.0*h)
		+ zs(i-1)*pow(xs(i)-x,3.0)/(6.0*h)
		+ (ys(i)/h-zs(i)*h/6.0)*(x-xs(i-1))
		+ (ys(i-1)/h-zs(i-1)*h/6.0)*(xs(i)-x);
	    return val;
	}
    }
    return ys.back(); // just to suppress warning
}

void check_range(double val, double val_min, double val_max,
		 string val_name, bool& gotError)
{
    if (val<val_min || val>val_max){
	cerr << "invalid " << val_name << " detected.\n";
	gotError = true;
    }
}

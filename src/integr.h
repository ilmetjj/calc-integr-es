#ifndef __INTEGR_H__
#define __INTEGR_H__

class integr;

#include <cmath>
#include <random>
#include <chrono>

using namespace std;

class integr
{
	double (*foo)(double);
	double a, b;
	int n;
	mt19937_64 eng;
public:
	integr(double _a, double _b, double (*_foo)(double));
	double rett(int n);
	double trap(int n);
	double montecarlo(int n, double ymax);
	double montecarlo(int n);
};

#endif
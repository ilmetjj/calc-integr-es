#ifndef __INTEGR_H__
#define __INTEGR_H__

class integr;

#include <cmath>
#include <random>
#include <chrono>
#include <gmp.h>

using namespace std;

void randgen(mpf_t& out, mpf_t a, mpf_t b, gmp_randstate_t state, mp_bitcnt_t prec);

class integr
{
	double a, b;
	double (*foo)(double);
	int n;
	mt19937_64 eng;
public:
	integr(double _a, double _b, double (*_foo)(double));
	double rett(int n);
	double trap(int n);
	double montecarlo(int n, double ymax);
	double montecarlo(int n);
};

class mp_integr
{
	mp_bitcnt_t prec;
	mpf_t a, b;
	void(*foo)(mpf_t& out, mpf_t in);
	gmp_randstate_t rndst;
public:
	mp_integr(mp_bitcnt_t _prec, mpf_t _a, mpf_t _b, void(*foo)(mpf_t& out, mpf_t in));
	void rett(mpf_t& out, mpz_t n);
	void trap(mpf_t& out, mpz_t n);
	void montecarlo(mpf_t& out, mpz_t n);
};

#endif
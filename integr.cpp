#include <iostream>
#include "integr.h"

using namespace std;

integr::integr(double _a, double _b, double (*_foo)(double)):a(_a),b(_b),foo(_foo){
	srand(time(NULL));
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	eng=mt19937_64(seed);
}
double integr::rett(int n){
	double S=0, l=(b-a)/n;

	for(int i=0; i<n; i++){
		S+=foo(a+i*l)*l;
	}

	return S;
}
double integr::trap(int n){
	double S=0, l=(b-a)/n;

	for(int i=0; i<n; i++){
		S+=(foo(a+i*l)+foo(a+i*l+l))*l/2;
	}

	return S;
}
double integr::montecarlo(int n){
	double S=0, l=(b-a)/n, t=0;
	for(int i=0; i<n; i++){
		
	}

	return S;
}


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
double integr::montecarlo(int n, double ymax){
	double S=0, A=2*(b-a)*ymax;
	int in=0, out=0;
	std::uniform_real_distribution<double> dis_x(a,b), dis_y(-ymax,ymax);
	double f, x, y;
	for(int i=0; i<n; i++){
		x=dis_x(eng);
		f=foo(x);
		y=dis_y(eng);
		if(abs(f)>ymax){
			cerr<<"f(x)=	"<<f<<"	>ymax=	"<<ymax<<endl;
		}
		if(f*y<0)
			out++;
		else{
			if(abs(f)<abs(y))
				out++;
			else
				in++;
		}
	}
	S=double(in)*double(A)/double(n);

	return S;
}
double integr::montecarlo(int n){
	double ymax=2*foo(a);

	for (double i = a+(b-a)/n; i <= b ; i+=(b-a)/n)
	{
		if(foo(i)>ymax)
			ymax=2*foo(i);
	}

	start:

	double S=0, A=2*(b-a)*ymax;
	int in=0, out=0;
	std::uniform_real_distribution<double> dis_x(a,b), dis_y(-ymax,ymax);
	double f, x, y;
	for(int i=0; i<n; i++){
		x=dis_x(eng);
		f=foo(x);
		y=dis_y(eng);
		if(abs(f)>ymax){
			cerr<<"f(x)=	"<<f<<"	>ymax=	"<<ymax<<endl;
			ymax=2*f;
			goto start;
		}
		if(f*y<0)
			out++;
		else{
			if(abs(f)<abs(y))
				out++;
			else
				in++;
		}
	}
	S=double(in)*double(A)/double(n);

	return S;
}


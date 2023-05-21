#include <iostream>
#include <fstream>
#include <iomanip>
#include "integr.h"

using namespace std;

double Sigmoid(double x) {
	if (x > -16 && x < 16)
		return exp(x) / (exp(x) + 1);
	else if (x >= 16)
		return 1;
	else
		return 0;
}

int main()
{

	cout << "Hello World!" << endl;
	integr s(0, 1, Sigmoid);

	ofstream fout;
	fout.open("plt");
	fout << "set logscale" << endl;
	fout << "plot 'inte.dat' u 1:2 w p, 'inte.dat' u 1:3 w p, 'inte.dat' u 1:4 w p" << endl;
	fout.close();

	fout.open("inte.dat");
	double re, tr, mc;
	for (long int i = 1; i < 1e10; i = i * 10) {
		re = s.rett(i);
		tr = s.trap(i);
		mc = s.montecarlo(i);
		fout << setprecision(1) << i << "	" << setprecision(64) << re << "	" << tr << "	" << mc << endl;
		cout << setprecision(1) << i << "	" << setprecision(64) << re << "	" << tr << "	" << mc << endl;

		//	if(i% int(1e6)==0)
		system("gnuplot plt -p");
	}
	fout.close();

	system("gnuplot plt -p");


	return 0;
}
#include <iostream>
#include <fstream>
#include <iomanip>
#include "integr.h"

using namespace std;

double Sigmoid(double x){
	if (x > -16 && x < 16)
		return exp(x) / (exp(x) + 1);
	else if (x >= 16)
		return 1;
	else
		return 0;
}

int main()
{
	cout<<"Hello World!"<<endl;
	integr s(0,1,Sigmoid);

	ofstream fout("inte.dat");
	for(int i=0; i<1e6; i+=1000){
		fout<<setprecision (1)<<i<<"	"<<setprecision (64)<<s.rett(i)<<"	"<<s.trap(i)<<endl;
		cout<<setprecision (1)<<i<<"	"<<setprecision (64)<<s.rett(i)<<"	"<<s.trap(i)<<endl;
	}
	fout.close();
	fout.open("plt");
	fout<<"set logscale"<<endl;
	fout<<"plot 'inte.dat' u 1:2 w d, 'inte.dat' u 1:3 w d"<<endl;
	fout.close();

	system("gnuplot plt -p");

	return 0;
}
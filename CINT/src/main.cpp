#include <iostream>
#include <fstream>
#include <iomanip>
#include "integr.h"
#include <gmpxx.h>
#include <thread>

using namespace std;

void foo(mpf_t& out, mpf_t x) {
	mpf_t temp1;
	mpf_init(temp1);

	mpf_pow_ui(temp1, x, 2);
	mpf_add_ui(temp1, temp1, 5);
	mpf_pow_ui(temp1, temp1, 3);
	mpf_sqrt(temp1, temp1);
	mpf_div(temp1, x, temp1);

	mpf_set(out, temp1);

	mpf_clear(temp1);
}

int main() {
	cout << "Function integrator" << endl;

	cout << "integrating the function 	f(x)=x/sqrt((x^2+5)^3)	 from -3 to 10" << endl;

	mp_bitcnt_t prec, lim;

	cout << "bit precision (es. 128):	"; cin >> prec;
	cout << "max 10 power (es. 9):	"; cin >> lim;

	mpf_set_default_prec(prec);
	int odpr = log10(pow(2, prec)) + 1;
	cout << "output digit precision =	~" << odpr << endl;
	cout.precision(odpr);


	mpf_t a, b;
	mpf_init_set_si(a, -3);
	mpf_init_set_si(b, 10);

	mp_integr I(prec, a, b, foo);

	mpf_t val;
	mpf_init(val);

	cout << "expected value (please fit with adequate precision):	"; mpf_inp_str(val, NULL, 10);
	mpf_out_str(NULL, 10, prec, val); cout << "	expected	" << val << endl;

	ofstream fout;
	fout.open("plt");
	fout << "set terminal png medium size 1920,1080" << endl;
	fout << "set output 'value.png'" << endl;
	fout << "set logscale" << endl;
	fout << "plot 'inte.dat' u 1:2 w p, 'inte.dat' u 1:3 w p, 'inte.dat' u 1:4 w p, 'inte.dat' u 1:5 w p" << endl;
	fout << "set output 'diffq.png'" << endl;
	fout << "set logscale" << endl;
	fout << "plot 'inte.dat' u 1:6 w p, 'inte.dat' u 1:7 w p, 'inte.dat' u 1:8 w p" << endl;
	fout << "set output 'time.png'" << endl;
	fout << "set logscale" << endl;
	fout << "plot 'inte.dat' u 1:9 w p, 'inte.dat' u 1:10 w p, 'inte.dat' u 1:11 w p" << endl;
	fout.close();

	mpf_t re, tr, mc;
	mpf_init(re);
	mpf_init(tr);
	mpf_init(mc);

	mpz_t _i, _max;
	mpz_init(_i);
	mpz_init_set_si(_max, 10);
	mpz_pow_ui(_max, _max, lim);

	fout.open("inte.dat");
	fout.precision(odpr);

	cout << "init for" << endl;

	for (mpz_set_ui(_i, 10); mpz_cmp(_i, _max) <= 0; mpz_mul_ui(_i, _i, 10)) {
		cout << _i << "	" << flush;
		fout << _i << "	" << flush;


		auto t1_start = std::chrono::high_resolution_clock::now();
		thread t1(&mp_integr::rett, &I, std::ref(re), _i);
		auto t2_start = std::chrono::high_resolution_clock::now();
		thread t2(&mp_integr::trap, &I, std::ref(tr), _i);
		auto t3_start = std::chrono::high_resolution_clock::now();
		thread t3(&mp_integr::montecarlo, &I, std::ref(mc), _i);

		t1.join();
		auto t1_end = std::chrono::high_resolution_clock::now();
		t2.join();
		auto t2_end = std::chrono::high_resolution_clock::now();
		t3.join();
		auto t3_end = std::chrono::high_resolution_clock::now();

		/*
				cerr << "ret ";
				I.rett(re, _i);
				cerr << " tr ";
				I.trap(tr, _i);
				cerr << " mc ";
				I.montecarlo(mc, _i);
				cerr << " end";
		*/
		cout << val << "	" << re << "	" << tr << "	" << mc << "	" << flush;
		fout << val << "	" << re << "	" << tr << "	" << mc << "	" << flush;

		mpf_sub(re, val, re);
		mpf_pow_ui(re, re, 2);
		mpf_sub(tr, val, tr);
		mpf_pow_ui(tr, tr, 2);
		mpf_sub(mc, val, mc);
		mpf_pow_ui(mc, mc, 2);

		cout << re << "	" << tr << "	" << mc << "	" << flush;
		fout << re << "	" << tr << "	" << mc << "	" << flush;

		cout << std::chrono::duration<double, std::milli>(t1_end - t1_start).count() << "	" << std::chrono::duration<double, std::milli>(t2_end - t2_start).count() << "	" << std::chrono::duration<double, std::milli>(t3_end - t3_start).count() << "	" << std::chrono::duration<double, std::milli>(t3_end - t1_start).count() << endl;
		fout << std::chrono::duration<double, std::milli>(t1_end - t1_start).count() << "	" << std::chrono::duration<double, std::milli>(t2_end - t2_start).count() << "	" << std::chrono::duration<double, std::milli>(t3_end - t3_start).count() << "	" << std::chrono::duration<double, std::milli>(t3_end - t1_start).count() << endl;

		cout << "Wall clock time passed: " << std::chrono::duration<double>(t3_end - t1_start).count() << " s\n" << endl;

		system("gnuplot plt -p");

	}
	fout.close();

	//	system("gnuplot plt -p");

	mpf_clear(a);
	mpf_clear(b);
	mpf_clear(val);
	mpf_clear(re);
	mpf_clear(tr);
	mpf_clear(mc);

	mpz_clear(_i);
	mpz_clear(_max);

	return 0;
}
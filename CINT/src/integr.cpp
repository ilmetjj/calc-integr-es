#include "integr.h"

using namespace std;

void randgen(mpf_t& out, mpf_t a, mpf_t b, gmp_randstate_t state, mp_bitcnt_t prec) {
	mpf_t temp1, temp2;
	mpf_init(temp1);
	mpf_init(temp2);

	mpf_urandomb(temp1, state, prec);
	mpf_sub(temp2, b, a);
	mpf_mul(temp1, temp1, temp2);
	mpf_add(temp1, temp1, a);

	mpf_set(out, temp1);

	mpf_clear(temp1);
	mpf_clear(temp2);
}

integr::integr(double _a, double _b, double (*_foo)(double)) :a(_a), b(_b), foo(_foo) {
	srand(time(NULL));
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	eng = mt19937_64(seed);
}
double integr::rett(int n) {
	double S = 0, l = (b - a) / n;

	for (int i = 0; i < n; i++) {
		S += foo(a + i * l) * l;
	}

	return S;
}
double integr::trap(int n) {
	double S = 0, l = (b - a) / n;

	for (int i = 0; i < n; i++) {
		S += (foo(a + i * l) + foo(a + i * l + l)) * l / 2;
	}

	return S;
}
double integr::montecarlo(int n, double ymax) {
	double S = 0, A = 2 * (b - a) * ymax;
	int in = 0, out = 0;
	std::uniform_real_distribution<double> dis_x(a, b), dis_y(-ymax, ymax);
	double f, x, y;
	for (int i = 0; i < n; i++) {
		x = dis_x(eng);
		f = foo(x);
		y = dis_y(eng);
		if (abs(f) > ymax) {
			cerr << "f(x)=	" << f << "	>ymax=	" << ymax << endl;
		}
		if (f * y < 0)
			out++;
		else {
			if (abs(f) < abs(y))
				out++;
			else
				in++;
		}
	}
	S = double(in) * double(A) / double(n);

	return S;
}
double integr::montecarlo(int n) {
	double ymax = 2 * abs(foo(a));

	for (double i = a + (b - a) / n; i <= b; i += (b - a) / n)
	{
		if (abs(foo(i)) > ymax)
			ymax = 2 * foo(i);
	}

start:

	double S = 0, A = 2 * (b - a) * ymax;
	int in = 0, out = 0;
	std::uniform_real_distribution<double> dis_x(a, b), dis_y(-ymax, ymax);
	double f, x, y;
	for (int i = 0; i < n; i++) {
		x = dis_x(eng);
		f = foo(x);
		y = dis_y(eng);
		if (abs(f) > ymax) {
			cerr << "f(x)=	" << f << "	>ymax=	" << ymax << endl;
			ymax = 2 * f;
			goto start;
		}
		if (f * y < 0)
			out++;
		else {
			if (abs(f) < abs(y))
				out++;
			else
				in++;
		}
	}
	S = double(in) * double(A) / double(n);

	return S;
}

mp_integr::mp_integr(mp_bitcnt_t _prec, mpf_t _a, mpf_t _b, void(*_foo)(mpf_t& out, mpf_t in)) : prec(_prec), foo(_foo) {
	mpf_set_default_prec(prec);

	mpf_init_set(a, _a);
	mpf_init_set(b, _b);

	gmp_randinit_mt(rndst);
	unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
	gmp_randseed_ui(rndst, seed);
}

void mp_integr::rett(mpf_t& out, mpz_t _n) {
	mpf_t S, l, n, i, temp1, temp2;
	mpz_t _i;

	mpf_init(S);
	mpf_init(l);
	mpf_init(n);
	mpf_init(i);
	mpf_init(temp1);
	mpf_init(temp2);

	mpz_init(_i);

	mpf_set_si(S, 0);
	mpf_set_si(l, 0);
	mpf_set_z(n, _n);

	mpf_sub(temp1, b, a);
	mpf_div(l, temp1, n);

	for (mpz_set_ui(_i, 0); mpz_cmp(_i, _n) < 0; mpz_add_ui(_i, _i, 1)) {
		mpf_set_z(i, _i);
		mpf_mul(temp1, i, l);
		mpf_add(temp1, temp1, a);
		foo(temp2, temp1);
		mpf_mul(temp1, temp2, l);
		mpf_add(S, S, temp1);
	}

	mpf_set(out, S);

	mpf_clear(S);
	mpf_clear(l);
	mpf_clear(n);
	mpf_clear(i);
	mpf_clear(temp1);
	mpf_clear(temp2);

	mpz_clear(_i);
}
void mp_integr::trap(mpf_t& out, mpz_t _n) {
	mpf_t S, l, n, i, temp1, temp2, temp3, temp4;
	mpz_t _i;

	mpf_init(S);
	mpf_init(l);
	mpf_init(n);
	mpf_init(i);
	mpf_init(temp1);
	mpf_init(temp2);
	mpf_init(temp3);
	mpf_init(temp4);

	mpz_init(_i);

	mpf_set_si(S, 0);
	mpf_set_si(l, 0);
	mpf_set_z(n, _n);

	mpf_sub(temp1, b, a);
	mpf_div(l, temp1, n);

	for (mpz_set_ui(_i, 0); mpz_cmp(_i, _n) < 0; mpz_add_ui(_i, _i, 1)) {
		mpf_set_z(i, _i);
		mpf_mul(temp1, i, l);
		mpf_add(temp1, temp1, a);
		foo(temp2, temp1);
		mpf_mul(temp3, i, l);
		mpf_add(temp3, temp3, a);
		mpf_add(temp3, temp3, l);
		foo(temp4, temp3);
		mpf_add(temp1, temp2, temp4);
		mpf_mul(temp1, temp1, l);
		mpf_div_ui(temp1, temp1, 2);
		mpf_add(S, S, temp1);
	}

	mpf_set(out, S);

	mpf_clear(S);
	mpf_clear(l);
	mpf_clear(n);
	mpf_clear(i);
	mpf_clear(temp1);
	mpf_clear(temp2);
	mpf_clear(temp3);
	mpf_clear(temp4);

	mpz_clear(_i);
}
void mp_integr::montecarlo(mpf_t& out, mpz_t _n) {
	mpf_t S, l, n, i, temp1, temp2, ymax, A, f, x, y;
	mpz_t _i, in, ot;

	mpf_init(S);
	mpf_init(l);
	mpf_init(n);
	mpf_init(i);
	mpf_init(temp1);
	mpf_init(temp2);
	mpf_init(ymax);
	mpf_init(A);
	mpf_init(f);
	mpf_init(x);
	mpf_init(y);

	mpz_init(_i);
	mpz_init(in);
	mpz_init(ot);

	mpf_set_z(n, _n);

	foo(ymax, a);
	mpf_abs(ymax, ymax);
	mpf_sub(temp1, b, a);
	mpf_div(l, temp1, n);
	mpf_add(temp1, a, l);

	for (mpf_set(i, temp1); mpf_cmp(i, b) <= 0; mpf_add(i, i, l))
	{
		foo(temp1, i);
		mpf_abs(temp1, temp1);
		if (mpf_cmp(temp1, ymax) > 0) {
			mpf_set(ymax, temp1);
		}
	}

start:

	mpf_set_ui(S, 0);

	mpf_sub(temp1, b, a);
	mpf_mul(temp1, temp1, ymax);
	mpf_mul_ui(A, temp1, 2);

	mpz_set_ui(in, 0);
	cerr << endl << in << endl;
	mpz_set_ui(ot, 0);
	cerr << endl << ot << endl;

	ofstream fout("point.dat");
	ofstream foutn("pointn.dat");
	ofstream foutb("pointo.dat");

	for (mpz_set_ui(_i, 0); mpz_cmp(_i, _n) < 0; mpz_add_ui(_i, _i, 1)) {
		randgen(x, a, b, rndst, prec);
		mpf_neg(temp1, ymax);
		randgen(y, temp1, ymax, rndst, prec);
		foo(f, x);

		mpf_abs(temp1, f);
		if (mpf_cmp(temp1, ymax) > 0) {
			cerr << "f(x)=	" << f << "	>ymax=	" << ymax << endl;
			mpf_set(ymax, temp1);
			goto start;
		}
		mpf_mul(temp2, f, y);
		if (mpf_cmp_si(temp2, 0) >= 0) {
			mpf_abs(temp2, y);
			if (mpf_cmp(temp1, temp2) >= 0) {
				if(mpf_cmp_si(y,0)>=0){
					mpz_add_ui(in, in, 1);
					//	cerr << endl << x << " " << y << " " << f << " " << in << endl;
					fout << x << "	" << y << "	" << f << endl;
				}
				else{
					mpz_sub_ui(in, in, 1);
					//	cerr << endl << x << " " << y << " " << f << " " << in << endl;
					foutn << x << "	" << y << "	" << f << endl;
				}
			}
			else
			{
				mpz_add_ui(ot, ot, 1);
				foutb << x << "	" << y << "	" << f << endl;
			}
		}
		else
		{
			mpz_add_ui(ot, ot, 1);
			foutb << x << "	" << y << "	" << f << endl;
		}


	}
	fout.close();
	foutb.close();
	fout.open("grp");
	fout << "set terminal png medium size 1920,1080" << endl;
	fout << "set output 'graph.png'" << endl;
	fout << "unset logscale" << endl;
	fout << "plot 'point.dat' u 1:2 w p, 'pointn.dat' u 1:2 w p, 'pointo.dat' u 1:2 w p, 'pointo.dat' u 1:3 w p" << endl;
	fout.close();

	system("gnuplot grp -p");

	cerr << "n	" << n << endl;
	cerr << "_n	" << _n << endl;
	cerr << "ot	" << ot << endl;
	cerr << "in	" << in << endl;
	mpf_set_z(temp1, in);
	cerr << "temp1	" << temp1 << endl;
	cerr << "A	" << A << endl;
	mpf_div(temp2, temp1, n);
	cerr << "temp2	" << temp2 << endl;
	mpf_mul(S, A, temp2);
	cerr << "S	" << S << endl;
	mpf_set(out, S);
	cerr << "out	" << out << endl;

	mpf_clear(S);
	mpf_clear(l);
	mpf_clear(n);
	mpf_clear(i);
	mpf_clear(temp1);
	mpf_clear(temp2);
	mpf_clear(ymax);
	mpf_clear(A);
	mpf_clear(f);
	mpf_clear(x);
	mpf_clear(y);

	mpz_clear(_i);
	mpz_clear(in);
}
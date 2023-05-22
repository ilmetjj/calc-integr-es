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
	cout << "max power of 10 exponent (es. 7):	"; cin >> lim;

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
	cout << "The expected value should be " << endl;
	cout << "0.169671234617571066749115017255743502594620151379863039692724881000094435110444653684750881767421905231456896782810353707446053147508095407713533504749271154765511607154456415080254568442457700513957746070996325209911170599369780380426016275370130025860989355890746917559816140369054180220025315311661909530382734861214246087094353142998705394623409713896335177722411254056264571013040273220949985780470011342995252997690491124994523292835336635737239634097221601121130163168152429527854133771190310521781601411922289978109543936481452119112326276290005205885108074912215744071507498856795834866527393778927029987223426728087430405513364694567485230270675786394370152052843240213309001840255009525455009056033208804050989517966011863157878929261591049844177267503600373749803847098294037829395378318179504758238521090332306587431356636966742536884533392366897393570458661714563633195957159998838164888966050778816498094124391246244960054149322076330194412450428999937060791148179075983086601983317521440867707610442849439454186475197601427106500886511792930086804237297376700261414104108242930290428476096819242233251975665156186857700410304795742159932856260798664008146039362539474978089212890925784287088887404781825025499802543532406759680609873496939873227451685966220083723654424941800162606456825794451603810839874292347192348665038736898020509493151590530554551236321686284052653230467932715501839409789140745962408080072994329858400062315168737360618255605839727986601422869292023948426074880653949801437474188609357649171776061681726124108289330762039928108857450101518134312169477939117672597556627231594819757829126736881190752319249606016938502192613788964837197834643434774319744697701234443713983814660916450285885749611127937193737868574128059724446762751771180620229116161737499009840598102836399332685292907442082634832762630185628742122202323821164318203996486760592504907256177821103831728296121475880307858334072" << endl;
	cout << "expected value (insert one or copy the above):	"; mpf_inp_str(val, NULL, 10);
	mpf_out_str(NULL, 10, prec, val); cout << "	expected	" << val << endl;

	ofstream fout;
	fout.open("plt");
	fout << "set terminal png medium size 1920,1080" << endl;

	fout << "set output 'value.png'" << endl;
	fout << "unset logscale" << endl;
	fout << "plot 'inte.dat' u 1:2 w l title 'value', 'inte.dat' u 1:3 w l title 'rett', 'inte.dat' u 1:4 w l title 'trap', 'inte.dat' u 1:5 w l title 'monc'" << endl;

	fout << "set output 'value-log.png'" << endl;
	fout << "set logscale" << endl;
	fout << "plot 'inte.dat' u 1:2 w p title 'value', 'inte.dat' u 1:3 w p title 'rett', 'inte.dat' u 1:4 w p title 'trap', 'inte.dat' u 1:5 w p title 'monc'" << endl;

	fout << "set output 'diffq.png'" << endl;
	fout << "unset logscale" << endl;
	fout << "plot 'inte.dat' u 1:6 w l title 'rett', 'inte.dat' u 1:7 w l title 'trap', 'inte.dat' u 1:8 w l title 'monc'" << endl;

	fout << "set output 'diffq-log.png'" << endl;
	fout << "set logscale" << endl;
	fout << "plot 'inte.dat' u 1:6 w p title 'rett', 'inte.dat' u 1:7 w p title 'trap', 'inte.dat' u 1:8 w p title 'monc'" << endl;

	fout << "set output 'time.png'" << endl;
	fout << "unset logscale" << endl;
	fout << "plot 'inte.dat' u 1:9 w l title 'rett', 'inte.dat' u 1:10 w l title 'trap', 'inte.dat' u 1:11 w l title 'monc'" << endl;

	fout << "set output 'time-log.png'" << endl;
	fout << "set logscale" << endl;
	fout << "plot 'inte.dat' u 1:9 w p title 'rett', 'inte.dat' u 1:10 w p title 'trap', 'inte.dat' u 1:11 w p title 'monc'" << endl;
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
		cout << _i << "	" << val << endl;
		fout << _i << "	" << val << flush;


		auto t1_start = std::chrono::high_resolution_clock::now();
		thread t1(&mp_integr::rett, &I, std::ref(re), _i);
		auto t2_start = std::chrono::high_resolution_clock::now();
		thread t2(&mp_integr::trap, &I, std::ref(tr), _i);
		auto t3_start = std::chrono::high_resolution_clock::now();
		thread t3(&mp_integr::montecarlo, &I, std::ref(mc), _i);

		t1.join();
		auto t1_end = std::chrono::high_resolution_clock::now();
		cout << " re= " << re << endl;
		t2.join();
		auto t2_end = std::chrono::high_resolution_clock::now();
		cout << " tr= " << tr << endl;
		t3.join();
		auto t3_end = std::chrono::high_resolution_clock::now();


		cout << " re= " << re << "	tr = " << tr << "	mc = " << mc << "	" << endl;
		fout << "	" << re << "	" << tr << "	" << mc << "	" << flush;

		mpf_sub(re, val, re);
		mpf_pow_ui(re, re, 2);
		mpf_sub(tr, val, tr);
		mpf_pow_ui(tr, tr, 2);
		mpf_sub(mc, val, mc);
		mpf_pow_ui(mc, mc, 2);

		cout << " re " << re << "	 tr " << tr << "	 mc " << mc << "	" << endl;
		fout << re << "	" << tr << "	" << mc << "	" << flush;

		cout << " re " << std::chrono::duration<double, std::milli>(t1_end - t1_start).count() << "	 tr " << std::chrono::duration<double, std::milli>(t2_end - t2_start).count() << "	 mc " << std::chrono::duration<double, std::milli>(t3_end - t3_start).count() << endl;
		fout << std::chrono::duration<double, std::milli>(t1_end - t1_start).count() << "	" << std::chrono::duration<double, std::milli>(t2_end - t2_start).count() << "	" << std::chrono::duration<double, std::milli>(t3_end - t3_start).count() << "	" << std::chrono::duration<double, std::milli>(t3_end - t1_start).count() << endl;

		cout << "Wall clock time passed: " << std::chrono::duration<double>(t3_end - t1_start).count() << " s\n" << endl;

		system("gnuplot plt -p");

	}
	fout.close();

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
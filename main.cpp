#include "FIR/FIRDesign.h"
#include "FIR/PolyphaseFilterbank.h"
#include "FIR/FractionalFIRFilter.h"
#include "NCO/NCOTable.h"
#include <iostream>

int main(){
	int rateI = 3;
	std::vector<double> taps{1,2,3,4,5,6};
	auto tapsf = Filter::FIR::Design::convertTapsForFilter<float>(taps, false, rateI);
	std::cout << tapsf->size() << std::endl;
	for(auto a: *tapsf){
		std::cout << a << std::endl;
	}

	return 0;
}

int main2(){
	int rateI = 5;
	auto taps = Filter::FIR::Design::buildFIRFilter(0.05, 0.1, 1, -140);
	//taps.push_back(0);
	//taps.push_back(0);
	//std::vector<double> taps{1,2,3,4,5,6};
	auto tapsf = Filter::FIR::Design::convertTapsForFilter<float>(taps, false, rateI);
	std::cout << tapsf->size() << std::endl;
	for(auto a: *tapsf){
	//	std::cout << a << std::endl;
	}

	Filter::FIR::FractionalFIR<float,std::complex<float>,std::complex<float>,false> fir(tapsf, rateI, 0.8);
	//Filter::FIR::FractionalFIR<float,float,float> fir(tapsf, rateI, 1.0);

	NCO::NCOTable<24, float> nco(1);

	//float f[50]={-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,1,1};
	//float ff[500];
	//fir.filter(f, 10, ff);
	//return 0;

	for(float f=0; f<0.5; f+=0.01){
		nco.setIncrement(nco.frequencyToIncrement(f));
		std::complex<float> out[10*2048];
		std::complex<float> buf[2048];
		nco.fillBuffer(buf, 2048);

		fir.filter(buf, 2048, out);

		float sum = 0;
		for (int i=1008; i<2008;i++){
			sum += std::abs(out[i]);
		}

		std::cerr << f<<","<<std::log(sum/1000.0*rateI)/std::log(10)*20<<","<<sum << std::endl;
	}
	return 0;
}

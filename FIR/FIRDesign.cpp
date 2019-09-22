#include <cmath>

#include "FIRDesign.h"
#include <iostream>

namespace Filter
{
namespace FIR
{
namespace Design
{

static double besselI0(double x)
{
	double result = 1;
	double fac = 1;
    x = x*x/4;
    double x2 = 1;

    for(int i=1; i<=15; i++) {
        fac *= i * i;
        x2 *= x;
        result += x2 / fac;
    }
    return result;
}

std::vector<double> buildKaiser(double rip, double tw, double fs)
{
    tw = 2*M_PI*tw/fs;

    if(rip > -21) {
        rip = -21;
    }

    unsigned int m = std::ceil((-7.95-rip)/(2.285*tw));

    double beta;
    if(rip < -50) {
        beta = 0.1102*(-8.7-rip);
    } else {
        beta = 0.07886*(-21.0-rip);
        beta += 0.5842*std::pow(-21.0-rip, 0.4);
    }


    double d = besselI0(beta);

    int windowLength = m+1;
    auto window = std::vector<double>(windowLength);

    for(unsigned int i=0; i<=m/2; i++) {
        double v = besselI0(beta*std::sqrt(1-std::pow(2*(double)i/(double)m-1,2)))/d;
        window[i] = v;
        window[m-i] = v;
    }

    return window;
}

std::vector<double> buildFIRBase(double fc, double fs, unsigned int taps, std::function<double(double)> gainFunc, double step)
{
    if(taps == 0) {
        return std::vector<double>(0);
    }

    auto result = std::vector<double>(taps);
    double ft = fc/fs;

    if(!gainFunc) {
        if(taps == 1) {
            result[0] = 1.0;
            return result;
        }

        /* If gain == 1 then we can solve the fourier integral analytically */
        for(unsigned int i=0; i<taps/2; i++) {
            double k = (double)i - (double)(taps-1)/2.0;
            double p = std::sin(2*M_PI*ft*k)/k/M_PI;
            result[i] = p;
            result[taps - i - 1] = p;
        }
        if(taps & 1) {
            result[taps/2] = 2*ft;
        }
    } else {
        if(taps == 1) {
            result[0] = gainFunc(0);
            return result;
        }

        step = (2*ft)/step;

        /* Calculate real discrete fourier transform */
        for(unsigned int i=0; i<=taps/2; i++) {
            double sum = gainFunc(0);
            for(double f = step; f <= ft; f+=step) {
                double gain = 2 * gainFunc(f*fs);
                sum += gain * std::cos(2*M_PI*f*(0.5+(double)i-(double)taps/2.0));
            }
            sum *= step;
            result[i] = sum;
            result[taps - i - 1] = sum;
        }
    }

    return result;
}

std::vector<double> buildFIRFilter(double fc, double tw, double fs, double rip,
                                  std::function<double(double)> gainFunc,
                                  double step)
{
    /* Build kaiser window that meets the specifications */
    auto window = buildKaiser(rip, tw, fs);
    auto filter = buildFIRBase(fc, fs, window.size(), gainFunc, step);

    /* Multiply window with filter */
    for(unsigned int i=0; i<filter.size(); i++) {
        filter[i] *= window[i];
    }

    return filter;
}

std::vector<double> buildFIRGaussian(unsigned int bits, unsigned int oversampling, double BT){
	unsigned int length = bits * oversampling + 1;

	std::vector<double> taps;

	double k = std::sqrt(2.0) * BT / std::sqrt(std::log(2.0));

	for(unsigned int i=0; i<length; i++){
		double t = (0.5 + (double)i - ((double)length)/2.0) / (double)oversampling;
		taps.push_back(std::exp(-t*t*M_PI*M_PI*k*k));
	}

	double div = 0;
	for(double tmp: taps){
		div += std::abs(tmp);
	}

	for(unsigned int i=0; i<length; i++){
		taps[i]/=div;
	}

	return taps;
}

}
}
}

#ifndef IIR_TRADITIONAL_CPP_
#define IIR_TRADITIONAL_CPP_

#include "IIRDesign.h"

using namespace std::complex_literals;


namespace Filter
{
namespace IIR
{
namespace Design
{

static std::vector<std::complex<double>> fillComplexVector(unsigned int len, double value = -1)
{
    std::vector<std::complex<double>> zeros(len);
    for(unsigned int i=0; i<zeros.size(); i++) {
        zeros[i] = value;
    }

    return zeros;
}

PoleZeroType butterworth(double fc, double fs, unsigned int order)
{
    std::vector<std::complex<double>> poles;

    double wc = Util::scaleFrequency(fc, fs);

    /* Calculate poles */
    for(unsigned int p=0; p<order/2 + order % 2; p++) {
        double theta = M_PI * (2*p + 1) / (2 * order);

        /* Calculate pole */
        std::complex<double> pole = wc * 1i * std::exp(theta * 1i);

        /* Transfor to z domain */
        pole = Util::transformBilinear(pole, fs);

        /* Store poles */
        poles.push_back(pole);

        if(p<order/2) {
            /* The final pole is real */
            poles.push_back(std::conj(pole));
        }
    }

    return std::make_tuple(poles, fillComplexVector(poles.size()), fc, fs);
}


PoleZeroType chebyshev(double fc, double fs, unsigned int order, double ripple)
{
    std::vector<std::complex<double>> poles;

    double wc = Util::scaleFrequency(fc, fs);

    /* Calculate poles */
    for(unsigned int p=0; p<order/2 + order % 2; p++) {
        double theta = M_PI * (2*p + 1) / (2 * order);

        /* Calculate pole at 1 rad/s */
        double epsilon = std::sqrt(std::pow(10, ripple/10) - 1);
        double ch = 1.0/order * std::asinh(1/epsilon);

        std::complex<double> pole = -std::sinh(ch) * std::sin(theta) * wc
                                    +std::cosh(ch) * std::cos(theta) * 1i * wc;

        pole = Util::transformBilinear(pole, fs);

        poles.push_back(pole);

        if(p<order/2) {
            /* The final pole is real */
            poles.push_back(std::conj(pole));
        }
    }

    return std::make_tuple(poles, fillComplexVector(poles.size()), fc, fs);
}

}
}
}


#endif /* IIR_TRADITIONAL_CPP_ */

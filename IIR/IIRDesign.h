#ifndef IIR_IIRDESIGN_H_
#define IIR_IIRDESIGN_H_

#include <tuple>
#include <vector>
#include <complex>

#include "IIRFilter.h"

namespace Filter
{
namespace IIR
{

using PoleZeroType = std::tuple<std::vector<std::complex<double>>,std::vector<std::complex<double>>, double, double>;

namespace Design
{

/* Functions to generate poles and zeroes */
PoleZeroType butterworth(double fc, double fs, unsigned int order);
PoleZeroType chebyshev(double fc, double fs, unsigned int order, double rippledB = 1);
PoleZeroType elliptic(double fc, double fs, unsigned int order, double passbandRippledB = 1, double stopbanddB = 80, double fstop = 0);

/* Function to make IIR filter */
std::vector<IIRCoefficients> buildIIR(const PoleZeroType& pz, bool highpass = false, bool useBiquads = false);

namespace Util
{
double scaleFrequency(double f, double fs);
std::complex<double> transformBilinear(std::complex<double> pole, double fs);
std::vector<double> expandTransformPoleZero(std::vector<std::complex<double>> poles, double transformCut = 0);
IIRCoefficients calculateDirect2(std::vector<std::complex<double>> cPoles, std::vector<std::complex<double>> cZeros, double transformCut = 0, double steadyGain = 1);
template <typename T> std::vector<T> polynomialMultiply(std::vector<T> a, std::vector<T> b)
{
    /*
     * a*x^2+bx+c
     * [c b a]
     */

    std::vector<T> result(a.size() + b.size() - 1);

    for(unsigned int i=0; i<a.size(); i++) {
        for(unsigned int j=0; j<b.size(); j++) {
            result[i + j] += a[i] * b[j];
        }
    }

    return result;
}

}
}
}
}



#endif /* IIR_IIRDESIGN_H_ */

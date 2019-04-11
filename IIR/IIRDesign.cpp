#include "IIRDesign.h"
#include <iostream>

namespace Filter
{
namespace IIR
{
namespace Design
{

double Util::scaleFrequency(double f, double fs)
{
    return 2.0 * M_PI * fs / M_PI * std::tan(M_PI * f / fs);
}

std::complex<double> Util::transformBilinear(std::complex<double> pole, double fs)
{
    /* Perform bilinear transform */
    return (1.0 + pole / (2 * fs)) / (1.0 - pole / (2 * fs));
}

std::vector<double> Util::expandTransformPoleZero(std::vector<std::complex<double>> poles, double transformCut)
{
    /*
     * Multiply together all poles as follows:
     * (z-p1)(z-p2)(...)(z-pn)
     */
    std::vector<std::complex<double>> combined = {1};

    for(unsigned int i=0; i<poles.size(); i++) {
        if(transformCut == 0) {
            combined = polynomialMultiply(combined, {-poles[i], 1});
        } else {
            double a = -std::cos(2*M_PI*transformCut);
            combined = polynomialMultiply(combined, {-a-poles[i], -1.0-poles[i]*a});
        }
    }

    std::vector<double> result (combined.size());
    for(unsigned int i=0; i<combined.size(); i++) {
        result[i] = std::real(combined[combined.size() - 1 - i]);
    }

    return result;
}


IIRCoefficients Util::calculateDirect2(std::vector<std::complex<double>> cPoles, std::vector<std::complex<double>> cZeros, double transformCut, double steadyGain)
{
    IIRCoefficients coeff;

    coeff.a = expandTransformPoleZero(cPoles, transformCut);
    coeff.b = expandTransformPoleZero(cZeros, transformCut);

    /* Calculate gain for z=-1 or z=1 */
    double sumZ = 0, sumP = 0, sign = 1;

    for(unsigned int i=0; i<coeff.a.size(); i++) {
        sumP += coeff.a[i] * sign;
        if(transformCut) {
            sign = -sign;
        }
    }

    sign = 1;

    for(unsigned int i=0; i<coeff.b.size(); i++) {
        sumZ += coeff.b[i] * sign;
        if(transformCut) {
            sign = -sign;
        }
    }

    steadyGain *= sumP / sumZ;

    coeff.finalGain = steadyGain;

    return coeff;
}

/*
 * This function requires that conjugate complex entries are provided together, so they get
 * packed into one biquad (otherwise the coefficients would be complex, and the rest of the code
 * does not handle that).
 *
 * For odd filters, the real pole goes last.
 */
std::vector<IIRCoefficients> buildIIR(const PoleZeroType& pz, bool highpass, bool useBiquads)
{
    std::vector<std::complex<double>> poles, zeros;
    double fc, fs;
    std::tie(poles, zeros, fc, fs) = pz;

    std::vector<IIRCoefficients> coeff;

    if(useBiquads) {
        for(unsigned int i=0; i<poles.size(); i+=2) {
            unsigned int remainingTerms = poles.size() - i;

            if(remainingTerms > 2) {
                remainingTerms = 2;
            }

            std::vector<std::complex<double>> subPoles(poles.begin() + i, poles.begin() + i + remainingTerms);
            std::vector<std::complex<double>> subZeros(zeros.begin() + i, zeros.begin() + i + remainingTerms);
            coeff.push_back(Util::calculateDirect2(subPoles, subZeros, highpass?(fc/fs):0));
        }
    } else {
        coeff.push_back(Util::calculateDirect2(poles, zeros, highpass?(fc/fs):0));
    }

    return coeff;
}

}
}
}

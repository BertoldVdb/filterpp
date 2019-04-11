#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <tuple>

#include "IIRDesign.h"

namespace Filter
{
namespace IIR
{
namespace Design
{

using namespace std::complex_literals;

static double agm(double a, double b, unsigned int iterations = 6)
{
    for(unsigned int i=0; i<iterations; i++) {
        double c = 0.5*(a+b);
        b = std::sqrt(a*b);
        a = c;
    }

    return a;
}

static double ellipK(double k)
{
    return M_PI/(2*agm(1, std::sqrt(1-k*k), 5));
}


static double jacobiTheta(unsigned int index, double v, double q, bool imaginary = false, unsigned int iterations = 10)
{
    /* Reference: http://mathworld.wolfram.com/JacobiThetaFunctions.html */
    /* Note: m == k^2 */
    /* Note: z = v*pi */
    /* If imaginary is set to true then v is interpret as j*v:
     * sin(x+j*y) = sin(x)*cosh(-y) - j*cos(x)*sinh(-y) =
     *              sin(x)*cosh(y)  + j*cos(x)*sinh(y) -=>
     * sin(j*y) = j*sinh(y)  */
    if(index == 1) {
        double result = 0;
        double sign = 1;

        for(unsigned int i=0; i<iterations; i++) {
            double a;
            if(!imaginary) {
                a = std::sin((2*i+1) * v * M_PI);
            } else {
                a = std::sinh((2*i+1) * v * M_PI); /* Warning, the result is also imaginary */
            }
            result += sign * std::pow(q, i*(i+1)) * a;
            sign = -sign;
        }

        return 2*std::pow(q, 1.0/4)*result;
    } else if(index == 4) {
        double result = 0;
        double sign = -1;

        for(unsigned int i=1; i<=iterations; i++) {
            double a;
            if(!imaginary) { /* The result is always real */
                a = std::cos(2*i * v * M_PI);
            } else {
                a = std::cosh(2*i * v * M_PI);
            }
            result += sign * std::pow(q, i*i) * a;
            sign = -sign;
        }

        return 1+2*result;
    }

    return 0;
}

static double calculateQ(double ka)
{
    /* Source: Elements of the Theory of Elliptic Functions, AMS, p115 */
    double l = (1-std::sqrt(ka))/(1+std::sqrt(ka))/2;
    return l + 2*std::pow(l,5) + 15*std::pow(l,9) + 150*std::pow(l,13);
}

static double minDegree(double k, double ka, double k1, double k1a)
{
    return ellipK(k) * ellipK(k1a) / ellipK(ka) / ellipK(k1);
}

/* This is the base formula from which ellipSNSpecific is derived to avoid useless
 * operations. If imaginary is true u is taken as u*j. The result is also imaginary.
 * It is kept here for documentation purposes.
static double ellipSNCanonical(double u, double m, bool imaginary = false, unsigned int iterations = 10){
	double k = std::sqrt(m);
	double ka = std::sqrt(1-k*k);

	double q = calculateQ(ka);
	double K = ellipK(k);
	double v = u/(2*K);

	return 1/std::sqrt(k) * jacobiTheta(1, v, q, imaginary, iterations) / jacobiTheta(4, v, q, imaginary, iterations);
}
*/

static double ellipSNSpecific(double v, double q, bool imaginary = false, unsigned int iterations = 10)
{
    return jacobiTheta(1, v, q, imaginary, iterations) / jacobiTheta(4, v, q, imaginary, iterations);
}

static double ellipZ(double k, double a)
{
    return std::sqrt((1+k*a*a)*(1+a*a/k));
}

static double findStopbandForOrder(unsigned int targetOrder, double k1, double k1a)
{
    /* Very simple bisection approach */
    double a = 1.00001, b = 1.99999, c, res;
    do {
        c = (a + b)/2;
        double k =  1/c;
        double ka = std::sqrt(1.0 - k*k);
        res = minDegree(k, ka, k1, k1a) - targetOrder;
        if(res > 0) {
            a=c;
        } else {
            b=c;
        }
    } while(std::abs(res) >= 1e-5 && std::abs(b-a) >= 1e-5);

    return b;
}

PoleZeroType elliptic(double fc, double fs, unsigned int order, double passbanddB, double stopbanddB, double fstop)
{
    double wc = Util::scaleFrequency(fc, fs);
    double ws;

    if(passbanddB > 5.5) {
        passbanddB = 5.5;
    }
    if(stopbanddB < 6.5) {
        stopbanddB = 6.5;
    }


    double k1 = std::sqrt((std::pow(10, passbanddB/10)-1)/(std::pow(10, stopbanddB/10)-1));
    double k1a = std::sqrt(1.0-k1*k1);

    if(order > 0) {
        ws = findStopbandForOrder(order, k1, k1a);
        ws *= wc;
    } else {
        ws = Util::scaleFrequency(fstop, fs);
    }

    double k =  wc/ws;
    double ka = std::sqrt(1.0-k*k);

    if(order == 0) {
        order = std::ceil(minDegree(k, ka, k1, k1a));
    }

    /* Calculate q once for all elliptic operations */
    double q = calculateQ(ka);

    double base = 0.5;
    if(order % 2) {
        /* Order is odd */
        base = 1;
    }

    /* Calculate p1 */
    double p = 1.0 / (2.0 * M_PI * order) * std::log(1.0+2.0/(std::pow(10, passbanddB/20) - 1));
    double p1 = ellipSNSpecific(p, q, true);

    /* Calculate poles and zeros */
    std::vector<std::complex<double>> poles;
    std::vector<std::complex<double>> zeros;
    double s = std::sqrt(k)/wc;
    for(unsigned int i=0; i<order/2; i++) {
        double index = base + i;
        double p2 = ellipSNSpecific(index/order, q);

        double z = 1.0/(s * p2);
        zeros.push_back(z * 1i);
        zeros.push_back(z * -1i);

        double real = -1*p1*ellipZ(-k, p2);
        double imag =    p2*ellipZ( k, p1);
        double scale = s*(1+p1*p1*p2*p2);
        poles.push_back((real+1i*imag)/scale);
        poles.push_back((real-1i*imag)/scale);
    }

    if(order % 2) {
        /* Odd filters have one more pole */
        poles.push_back(-p1 / s);
    }

    for(auto& pole: poles) {
        pole = Util::transformBilinear(pole, fs);
    }
    for(auto& zero: zeros) {
        zero = Util::transformBilinear(zero, fs);
    }

    /* Equal amount of poles and zeros */
    while(zeros.size() < poles.size()) {
        zeros.push_back(-1);
    }

    return std::make_tuple(poles, zeros, fc, fs);
}

}
}
}

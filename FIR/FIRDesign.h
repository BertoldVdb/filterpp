#ifndef FIR_FIRDESIGN_H_
#define FIR_FIRDESIGN_H_

#include <tuple>
#include <functional>
#include <vector>
#include <memory>

namespace Filter
{
namespace FIR
{
namespace Design
{

std::vector<double> buildKaiser(double rip, double tw, double fs);
std::vector<double> buildFIRBase(double fc, double fs, unsigned int taps, std::function<double(double)> gainFunc = nullptr, double step = 200);

std::vector<double> buildFIRFilter(double fc,
									double tw = 0.1f,
									double fs = 1.0f,
									double rip = -60,
									std::function<double(double)> gainFunc = nullptr,
									double step = 200);

std::vector<double> buildFIRGaussian(unsigned int bits, unsigned int oversampling, double BT);

template <typename T> std::shared_ptr<std::vector<T>> convertTapsForFilter(std::vector<double> in, bool padFront, unsigned int rateI = 1, double gain = 1.0){
	/* Pad taps with zeros */
	unsigned int i;
	for(i=0; i < in.size(); i+=rateI){}
	in.resize(i);

	/* Also add rateI zeros upfront to we can interpolate into the 'future' */
	auto out = std::make_shared<std::vector<T>>(((padFront?rateI:0) + in.size()) * 2);

	auto delayLen = out->size()/2/rateI;
	for(unsigned int j=0; j<rateI; j++){
		for(unsigned int i=0; i<delayLen; i++){
			T value = in[i*rateI+j];
			if(padFront){
				value = (!i)?0:in[(i-1)*rateI+j];
			}
			(*out)[i+j*delayLen*2] = (*out)[i+j*delayLen*2+delayLen] = value;
		}
	}

	return out;
}

}
}
}


#endif /* FIR_FIRDESIGN_H_ */

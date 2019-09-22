#ifndef FILTER_CIC_CIC_H_
#define FILTER_CIC_CIC_H_

#include <iostream>
#include "IntComb.h"
#include <array>
#include <cmath>

namespace Filter
{
namespace CIC
{

template <typename inType, typename outType, typename internalType = int32_t, unsigned int maxOrder = 16> class CIC
{
public:
    CIC() {}

    CIC(unsigned int order, unsigned int rate, internalType internalGain = 1):
        order_(order),
        rate_(rate),
		internalGain_(internalGain)
    {
        gain_ = 1;
        for(unsigned int i=0; i<order; i++) {
            gain_ *= rate;
        }
    }

    internalType getGrowth() const
    {
        return gain_;
    }

    double getGain(double outputFrequency)
    {
        if(outputFrequency == 0) {
            return 1;
        }

        double tmp = std::sin(M_PI * (double)rate_ * outputFrequency);
        tmp /= std::sin(M_PI * outputFrequency);
        tmp /= rate_;
        tmp = std::abs(tmp);
        auto a = std::pow(tmp, order_);

        return a;
    }

    size_t filter(inType* samplesIn, size_t numSamples, outType* samplesOut, unsigned int incIn = 1, unsigned int incOut = 1)
    {
    	size_t numOut = 0;
    	size_t cntIn = 0;
    	size_t cntOut = 0;

        for(size_t i=0; i<numSamples; i++) {
            internalType sample = samplesIn[cntIn] * internalGain_;
            cntIn+=incIn;

            /* Run integrator stages */
            for(unsigned int j=0; j<order_; j++) {
                sample = int_[j].update(sample);
            }

            outIndex_++;
            if(outIndex_ >= rate_) {
                outIndex_ = 0;

                /* Run comb stages */
                for(unsigned int j=0; j<order_; j++) {
                    sample = comb_[j].update(sample);
                }

                samplesOut[cntOut] = sample;
                cntOut+=incOut;
                numOut++;
            }
        }
        return numOut;
    }

    void reset(){
    	for(auto& i: int_){
    		i.reset();
    	}
    	for(auto& i: comb_){
			i.reset();
		}

    	outIndex_ = 0;
    }

private:
    unsigned int order_ = 0, rate_ = 1, outIndex_ = 0;

    internalType gain_;
    inType internalGain_;

    std::array<Integrator<internalType>, maxOrder> int_;
    std::array<Comb<internalType>, maxOrder> comb_;
};


}
}

#endif /* FILTER_CIC_CIC_H_ */

#include <cmath>
#include <cstdint>

#include "../CIC/CIC.h"
#include "../FIR/FIRDesign.h"
#include "../FIR/FractionalFIRFilter.h"
#include "../Util/Math.h"

#include "Resampler.h"

#include <iostream>
#include <functional>
#include <stdexcept>
#include <memory>

#ifndef FILTER_RESAMPLER_DOWNSAMPLER_H_
#define FILTER_RESAMPLER_DOWNSAMPLER_H_

namespace Filter
{
namespace Resampler
{

template <typename inT, typename outT, typename cicInternalType = int32_t, typename internalT = float> class Downsampler
{
public:
    Downsampler(ResamplerConfig<internalT>& cfg):
        cfg_(cfg)
    {
        /*
        * The CIC filter makes half band filters, so we use it
        * to get the rate to 2-4x of the output rate
        */
        unsigned int cicD = 1;
        while((cfg_.inputRate / cicD) >= cfg_.cicOversampleFactor*cfg_.outputRate) {
            cicR_++;
            cicD *= 2;
        }

        if(cicR_ == 1) {
            /* Avoid useless work */
            cfg_.cicOrder = 0;
        }

        double firInputRate = cfg.inputRate / cicD;

        /* Create CIC core */
        cic_ = CIC::CIC<inT, internalT, cicInternalType>(cfg_.cicOrder, cicD, cfg_.cicBoost);
        std::function<double(double)> cicGainFunc = [this](double f) {
        	/* Convert the frequency to a fraction of the input rate */
        	f /= cfg_.inputRate;

        	return 1.0/cic_.getGain(f)/cic_.getGrowth()/cfg_.cicBoost;
        };
        actualOutputRate_ = cfg_.inputRate / (float)cicD;

        if(cicR_ == 1) {
            /* If we don't supply the gain function the integral is solved analytically */
            cicGainFunc = nullptr;
        }

        unsigned long interpolation;
        double decimation;

        if(cfg_.rational) {
            interpolation = cfg_.outputRate;
            decimation = firInputRate;

            /* Reduce fraction */
            unsigned long d = Util::gcd((unsigned long)decimation, interpolation);
            decimation /= d;
            interpolation /= d;

            if(interpolation > cfg_.maxInterpolation) {
                throw std::runtime_error("Interpolation factor exceeds maximum!");
            }

            actualOutputRate_*= interpolation;
            actualOutputRate_/= decimation;

        }else{
            interpolation = cfg_.maxInterpolation;
            decimation = firInputRate * interpolation / cfg_.outputRate;

            actualOutputRate_*= interpolation;
            actualOutputRate_/= decimation;

        }

		if(!cfg.taps){
			/* Create FIR core */
			auto taps = FIR::Design::buildFIRFilter(cfg_.cutoffFrequency,
													cfg_.transitionWidth,
													firInputRate * interpolation,
													cfg_.desiredRipple,
													cicGainFunc);

			for(auto& tap: taps) {
				tap *= interpolation;
			}
			cfg.taps = FIR::Design::convertTapsForFilter<internalT>(taps, true, interpolation);
		}

		fractionalFir_ = FIR::FractionalFIR<internalT, internalT, outT>(cfg.taps, cfg.filterInterpolate, interpolation, decimation);
    }

    size_t filter(inT* samplesIn, size_t numSamples, outT* samplesOut, unsigned int incIn = 1, unsigned int incOut = 1)
    {
        if(cicR_ != 1) {
        	internalT samplesCIC[numSamples/cicR_+1];
            auto samplesAfterCIC = cic_.filter(samplesIn, numSamples, samplesCIC, incIn, 1);

            return fractionalFir_.filter(samplesCIC, samplesAfterCIC, samplesOut, 1, incOut);
        } else {
            return fractionalFir_.filter(samplesIn, numSamples, samplesOut, incIn, incOut);
        }
    }

    double getActualOutputRate()
    {
        return actualOutputRate_;
    }

    void reset(){
    	cic_.reset();
    	fractionalFir_.reset();
    }

private:
    CIC::CIC<inT, internalT, cicInternalType> cic_;
    FIR::FractionalFIR<internalT, internalT, outT> fractionalFir_;

    ResamplerConfig<internalT> cfg_;
    double actualOutputRate_;

    unsigned int cicR_ = 1;
};

}
} /* namespace Filter */

#endif /* FILTER_RESAMPLER_DOWNSAMPLER_H_ */

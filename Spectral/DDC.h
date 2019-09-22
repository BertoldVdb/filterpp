#ifndef FILTERPP_UTIL_DDC_H_
#define FILTERPP_UTIL_DDC_H_

#include <complex>
#include "../Resampler/Downsampler.h"
#include "../NCO/NCOTable.h"

namespace Spectral {

template <int ncoPhaseBits = 12, typename ncoPhaseAccType = uint32_t, typename cicInternalType = uint32_t, typename resamplerInternalType = float> class DDC {
public:
	DDC(Filter::Resampler::ResamplerConfig<resamplerInternalType>& config):
				resamplerI_(config),
				resamplerQ_(config),
				nco_(config.inputRate, 0){


	}

	void setReceiveFrequency(float freq){
		nco_.setIncrement(nco_.frequencyToIncrement(freq));
	}

	inline size_t processInput(std::complex<float>* in, size_t numIn, std::complex<float>* out){
		if(!numIn){
			return 0;
		}

		std::complex<float> bufferNco[numIn];
		nco_.fillBuffer(bufferNco, numIn);

		/* Mix input signal with NCO */
		for(size_t i=0; i<numIn; i++){
			bufferNco[i] *= in[i];
		}

		lastNCO_ = bufferNco[numIn-1];

		/* Run resampling filter */
		return resampleBuffer(bufferNco, numIn, out);
	}

	std::complex<float> getLastNCO(){
		return lastNCO_;
	}

	void reset(){
		nco_.clearAcc();
		resamplerI_.reset();
		resamplerQ_.reset();
	}

private:
	inline size_t resampleBuffer(std::complex<float>* in, size_t numIn, std::complex<float>* out){
		resamplerI_.filter((float*)in, numIn, (float*)out, 2, 2);
		return resamplerQ_.filter(((float*)in) + 1, numIn, ((float*)out) + 1, 2, 2);
	}

	Filter::Resampler::Downsampler<float, float, cicInternalType> resamplerI_;
	Filter::Resampler::Downsampler<float, float, cicInternalType> resamplerQ_;

	NCO::NCOTable<ncoPhaseBits, float, ncoPhaseAccType> nco_;

	std::complex<float> lastNCO_;
};

}

#endif /* FILTERPP_UTIL_DDC_H_ */

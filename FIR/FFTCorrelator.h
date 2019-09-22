#ifndef FFTCORRELATOR_H_
#define FFTCORRELATOR_H_

#include <ffts/ffts.h>
#include <complex>
#include <vector>
#include "Keeper.h"
#include <iostream>

namespace Filter
{
namespace FIR
{

template <typename T> class FFTCorrelator {
public:
	size_t processSamples(std::complex<float>* reference, int shift, T* out, T* in, size_t numIn){
		return processSamplesMultiple(&reference, &shift, 1, out, in, numIn);
	}

	void clear(){
		correlationBuffer_.clear();
	}

	FFTCorrelator(unsigned int blockSize, unsigned int numCoeff):
			blockSize_(blockSize),
			numCoeff_(numCoeff),
			correlationBuffer_(blockSize_){

		if (std::is_same<T, std::complex<float>>::value){
			fftPlanFoward_ = ffts_init_1d(blockSize, FFTS_FORWARD);
			fftPlanBackward_ = ffts_init_1d(blockSize, FFTS_BACKWARD);
		}else{
			fftPlanFoward_ = ffts_init_1d_real(blockSize, FFTS_FORWARD);
			fftPlanBackward_ = ffts_init_1d_real(blockSize, FFTS_BACKWARD);
		}
	}

	~FFTCorrelator() {
		if(fftPlanFoward_){
			ffts_free(fftPlanFoward_);
		}
		if(fftPlanBackward_){
			ffts_free(fftPlanBackward_);
		}
	}

	std::complex<float>* newReferenceVector(T* coeff){
		/* Zero pad reference */
		T fftIn[blockSize_] = {};
		std::memcpy(fftIn, coeff, numCoeff_*sizeof(T));

		/* Calculate FFT */
		auto fftOutput = new std::complex<float>[blockSize_];
		ffts_execute(fftPlanFoward_, fftIn, fftOutput);

		/* Prepare points */
		for(unsigned int i = 0; i<blockSize_; i++){
			fftOutput[i] = std::conj(fftOutput[i])/(float)blockSize_;
		}

		return fftOutput;
	}

	size_t processSamplesMultiple(std::complex<float>** reference, int* shift, unsigned numReference, T* out, T* in, size_t numIn){
		size_t offset = 0;

		size_t samplesWritten = 0;
		auto step = blockSize_ - numCoeff_ + 1;

		while(numIn){
			auto rem = correlationBuffer_.remaining();

			if(in){
				if(rem > numIn){
					rem = numIn;
				}
				correlationBuffer_.insert(in + offset, rem);
			}else{
				T zeros[rem]={};
				correlationBuffer_.insert(zeros, rem);
				numIn = 0;
				rem = 0;
			}

			if(correlationBuffer_.remaining() == 0){
				T fftBuf0[blockSize_];
				correlationBuffer_.extract(fftBuf0);

				std::complex<float> fftBuf1[blockSize_];
				ffts_execute(fftPlanFoward_, fftBuf0, fftBuf1);

				for(unsigned int j=0; j<numReference; j++){
					auto ref = reference[j];
					auto s = shift[j];

					int important = (std::is_same<T, std::complex<float>>::value)?blockSize_:(blockSize_/2+1);

					if(!s){
						for(int i=0; i<important; i++){
							fftBuf1[i] *= ref[i];
						}
					}else{
						for(int i=0; i<important; i++){
							int refIndex = i + s;
							if(refIndex < 0){
								refIndex += blockSize_;
							}
							if((unsigned int)refIndex >= blockSize_){
								refIndex -= blockSize_;
							}

							if(std::is_same<T, std::complex<float>>::value){
								fftBuf1[i] *= ref[refIndex];
							}else{
								if(refIndex >= important){
									fftBuf1[i] *= std::conj(ref[blockSize_-refIndex]);
								}else{
									fftBuf1[i] *= ref[refIndex];
								}
							}
						}
					}

					ffts_execute(fftPlanBackward_, fftBuf1, fftBuf0);

					std::memcpy(out+samplesWritten, fftBuf0, step*sizeof(T));
					samplesWritten += step;
				}

				correlationBuffer_.drop(step);
			}

			offset += rem;
			numIn -= rem;
		}

		return samplesWritten;
	}

private:
	unsigned int blockSize_;
	unsigned int numCoeff_;

	Util::Keeper<T> correlationBuffer_;

	ffts_plan_t* fftPlanFoward_ = nullptr;
	ffts_plan_t* fftPlanBackward_ = nullptr;

};

}
}

#endif /* FFTCORRELATOR_H_ */

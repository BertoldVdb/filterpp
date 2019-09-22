#ifndef FILTER_FIR_FRACTIONALFIRFILTER_H_
#define FILTER_FIR_FRACTIONALFIRFILTER_H_

#include <vector>
#include <iostream>

namespace Filter
{
namespace FIR
{

template <typename tapType, typename inType, typename outType> class FractionalFIR
{
public:
	FractionalFIR(){

	}

    FractionalFIR(std::shared_ptr<std::vector<tapType>> taps, bool interpolate = true, unsigned int rateI = 1, float rateD = 1.0f):
        taps_(taps),
		interpolate_(interpolate),
        rateI_(rateI),
        rateD_(rateD)
    {

    	/* Check tabs */
    	if(taps->size() % rateI){
			throw std::runtime_error("Taps not correctly padded");
		}

        reset();
    }

    void setDecimationFactor(float rateD)
    {
        rateD_ = rateD;
    }

    /*
     * This function allows resetting the filter phase, reducing computational load
     * if decimation is integer.
     */
    void resetPhase(){
    	outIndex_ = std::floor(outIndex_);
    }

    size_t filter(inType* samplesIn, size_t numSamples, outType* samplesOut, unsigned int incIn = 1, unsigned int incOut = 1)
    {
    	size_t numOut = 0;
    	size_t cntIn = 0;
        size_t cntOut = 0;

        if(!taps_){
        	return 0;
        }

		for(unsigned int i=0; i<numSamples; i++) {
			delayLine_[delayIndex_] = samplesIn[cntIn];
			cntIn += incIn;

			while(outIndex_ < rateI_) {
				/* Interpolate between two polyphase taps */
				float o = calcPolyphase(0);

		        bool doInterpolate = interpolate_ && outIndex_ != std::floor(outIndex_);

				if(doInterpolate){
					//TODO: Check if this is fully correct, I don't use it...
					float o2 = calcPolyphase(1);

					float remainder = outIndex_ - std::floor(outIndex_);
					o = remainder * o2 + (1-remainder) * o;
				}

				outIndex_ = outIndex_ + rateD_;
				samplesOut[cntOut] = o;
				cntOut += incOut;

				numOut++;
			}

			outIndex_ = outIndex_ - rateI_;

			if(delayIndex_){
				delayIndex_ --;
			}else{
				delayIndex_ = delayLine_.size() - 1;
			}
		}

        return numOut;
    }

    void initDelayLine(inType value){
     	for(unsigned int j=0; j<delayLine_.size(); j++) {
     		delayLine_[j] = value;
     	}
    }

    void reset(){
		outIndex_ = 0;
		delayLine_ = std::vector<inType>(taps_->size() / rateI_ / 2);
		delayIndex_ = delayLine_.size() - 1;
    }

private:
    inline outType calcPolyphase(unsigned int offset)
    {
        unsigned int fIndex = delayIndex_ + 1;
        unsigned int outIndex = outIndex_;

        outIndex += offset;
        if(outIndex >= rateI_){
        	outIndex -= rateI_;

        	/* Use the next sample */
        	fIndex--;
        }

        auto end = delayLine_.size();
        auto tapOffset = end - fIndex + 2*outIndex*end;

		auto taps = taps_->data() + tapOffset;
		auto delayLine = delayLine_.data();
		auto o = MAC(0, delayLine, taps, end);

        return o;
    }

    inline outType MAC(outType o, inType* a, tapType* b, unsigned int end){
		for(unsigned int j=0; j<end; j++) {
			o = o + a[j]*b[j];
		}
    	return o;
    }

    std::shared_ptr<std::vector<tapType>> taps_;
    bool interpolate_ = true;
    unsigned int rateI_ = 1;
    float rateD_ = 1;
    std::vector<inType> delayLine_;

    unsigned int delayIndex_ = 0;
    float outIndex_ = 0;
};

} /* namespace Fir */
} /* namespace Filter */



#endif /* FILTER_FIR_FRACTIONALFIRFILTER_H_ */

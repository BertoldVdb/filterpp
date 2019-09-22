#include "GFSKModulator.h"

#include <iostream>
#include <cstring>

namespace Mod {
namespace GFSK{

GFSKModulator::GFSKModulator(unsigned int filterBits, float BT, unsigned int oversampling, float modIndex):
			filterBits_(filterBits),
			BT_(BT),
			oversampling_(oversampling),
			modIndex_(modIndex),
			gfskFilter_(Filter::FIR::Design::convertTapsForFilter<float>(Filter::FIR::Design::buildFIRGaussian(filterBits, oversampling, BT), false)){

}

size_t GFSKModulator::calculateInstantenousFrequency(float* output, bool* bits, size_t numBits, size_t encodeLen, bool padRepeat){
	auto len = (numBits + filterBits_) * oversampling_;
	float modulatorInput[len]={};

	/* Calculate fixed scale factor */
	float fixed = (float)M_PI*modIndex_/(float)oversampling_;

	/* Load bits into filter */
	float last = 0;
	for(size_t i=0; i<numBits; i++){
		for(unsigned j=0; j<oversampling_; j++){
			last = modulatorInput[oversampling_*i+j] = fixed * (bits[i]?1.0f:-1.0f);
		}
	}

	if(padRepeat){
		/* Fill filter with first sample */
		gfskFilter_.initDelayLine(modulatorInput[0]);
	}else{
		last = 0;
		gfskFilter_.initDelayLine(0);
	}

	/* Repeat last sample */
	for(size_t i=numBits*oversampling_; i<(numBits+filterBits_)*oversampling_; i++){
		modulatorInput[i] = last;
	}

	/* Calculate result */
	gfskFilter_.filter(modulatorInput, len, modulatorInput);

	auto loss = (filterBits_*oversampling_)/2;

	/* If length is zero, use the maximum */
	if(!encodeLen){
		encodeLen = -1;
	}
	if(encodeLen > len-2*loss){
		encodeLen = len-2*loss;
	}

	auto encodeStart = len-2*loss - encodeLen;
	auto offset = loss+encodeStart;

	/* Copy result to output */
	std::memcpy(output, modulatorInput+offset, encodeLen * sizeof(float));

	return encodeLen;
}

}
}

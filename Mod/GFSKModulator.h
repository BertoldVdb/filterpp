#ifndef GFSKMODULATOR_H_
#define GFSKMODULATOR_H_

#include <complex>
#include "../FIR/FIRDesign.h"
#include "../FIR/FractionalFIRFilter.h"


namespace Mod {
namespace GFSK{


class GFSKModulator {
public:
	GFSKModulator(unsigned int filterBits, float BT=0.5, unsigned int oversampling = 8, float modIndex = 0.5);
	~GFSKModulator();

	size_t calculateInstantenousFrequency(float* output, bool* bits, size_t numBits, size_t encodeLen = 0, bool padRepeat = true);

private:
	unsigned int filterBits_;
	float BT_;
	unsigned int oversampling_;
	float modIndex_;

	Filter::FIR::FractionalFIR<float, float, float> gfskFilter_;
};

}
}

#endif /* GFSKMODULATOR_H_ */

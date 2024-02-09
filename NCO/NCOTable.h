#ifndef NCOTABLE_H_
#define NCOTABLE_H_

#include <cmath>
#include <cstdint>
#include <random>
#include <atomic>
#include <thread>
#include <complex>

namespace NCO {

template <int phaseBits = 8, typename outType = float, typename phaseAccType = uint32_t, int gain = 1> class NCOTable {
public:
	NCOTable(float sampleRate = 1, int ditherShift = -1):
			sampleRate_(sampleRate),
			ditherShift_(ditherShift){

		do{
			dither_ = random();
		}while(dither_ == 0);

		initSineRom();
	}

	inline void fillBuffer(outType* sine, outType* cosine, unsigned int numSamples, unsigned int inc = 1){
		for(unsigned int i=0; i<numSamples; i++){
			phaseAccType phase = acc_ + phase_;

			if(ditherShift_ >= 0){
				dither_ = dither(dither_);
				phase += dither_ >> (phaseBits - ditherShift_);
			}

			phase >>= (sizeof(phaseAccType) * 8) - phaseBits;

			if(sine){
				*sine = getValueSin(phase);
				sine+=inc;
			}
			if(cosine){
				*cosine = getValueCos(phase);
				cosine+=inc;
			}

			acc_ += increment_;
		}
	}

	inline void fillBuffer(std::complex<outType>* out, unsigned int numSamples){
		fillBuffer((outType*)out+1, ((outType*)out), numSamples, 2);
	}

	inline void setIncrement(phaseAccType increment){
		increment_ = increment;
	}

	inline void clearAcc(){
		acc_ = 0;
	}

	inline void setPhase(phaseAccType phase){
		phase_ = phase;
	}

	inline phaseAccType frequencyToIncrement(float freq){
		float t = freq/sampleRate_;
		return t * (1ULL<<(sizeof(phaseAccType)*8));
	}

private:
	float sampleRate_;

	phaseAccType acc_ = 0;
	phaseAccType phase_ = 0;
	phaseAccType increment_ = 0;

	static const int sineRomSize = (1ULL<<(phaseBits-2)) + 1;

	void initSineRom(){
		if(sineRomDone_.load(std::memory_order_acquire)) {
			return;
		}

		/* Check if we need to make it */
		if(!sineRomFirst_.test_and_set()) {
			for(int i=0; i<sineRomSize; i++){
				sineRom_[i] = (float)gain * std::sin(M_PI/2.0/(float)(sineRomSize-1) * (float)i);
			}

			sineRomDone_.store(true, std::memory_order_release);
		}else{
			while(!sineRomDone_.load(std::memory_order_acquire)) {
				std::this_thread::yield();
			}
		}
	}

	inline uint32_t dither(uint32_t x){
		x ^= x << 13;
		x ^= x >> 17;
		x ^= x << 5;
		return x;
	}

	inline outType getValueSin(phaseAccType phase){
		return getRom(phase);
		//return gain*std::sin((float)phase/65536.0f*2.0*M_PI);
	}

	inline outType getValueCos(phaseAccType phase){
		return getRom(phase + sineRomSize - 1);
		//return gain*std::cos((float)phase/65536.0f*2.0*M_PI);
	}

	inline outType getRom(phaseAccType phase){
		auto romIndex = phase & (sineRomSize - 2);
		auto quadrant = (phase >> (phaseBits - 2)) & 3;

		switch(quadrant){
		case 0:
			return sineRom_[romIndex];
		case 1:
			return sineRom_[sineRomSize-1-romIndex];
		case 2:
			return -sineRom_[romIndex];
		default:
			return -sineRom_[sineRomSize-1-romIndex];
		}
	}

	int ditherShift_;
	uint32_t dither_;

	static outType sineRom_[sineRomSize ];
    static std::atomic_flag sineRomFirst_;
    static std::atomic<bool> sineRomDone_;
};

template <int phaseBits, typename outType, typename phaseAccType, int gain> std::atomic<bool> NCOTable<phaseBits, outType, phaseAccType, gain>::sineRomDone_ (false);
template <int phaseBits, typename outType, typename phaseAccType, int gain> std::atomic_flag NCOTable<phaseBits, outType, phaseAccType, gain>::sineRomFirst_;
template <int phaseBits, typename outType, typename phaseAccType, int gain> outType NCOTable<phaseBits, outType, phaseAccType, gain>::sineRom_[sineRomSize ];

}

#endif /* NCOTABLE_H_ */

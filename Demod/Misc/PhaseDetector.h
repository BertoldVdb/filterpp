#ifndef PHASEDETECTOR_H_
#define PHASEDETECTOR_H_

namespace Demod {

namespace Misc {

template <typename InType = float, typename WorkType = float> class ComplexPhaseDetector{
public:
	ComplexPhaseDetector(){
	}

	inline void processBuffer(InType* I1, InType* Q1, InType* I2, InType* Q2, InType* out, unsigned int numSamples){
		for(unsigned int k=0; k<numSamples; k++){
			WorkType i1 = I1[k];
			WorkType q1 = Q1[k];
			WorkType i2 = I2[k];
			WorkType q2 = Q2[k];

			//TODO: Scale by signal power, make >>15 depend on used types
			*out = (i1*i2+q1*q2) >> 15;
			out++;
		}
	}
private:
};

}

}

#endif /* PHASEDETECTOR_H_ */

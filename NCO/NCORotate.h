/*
 * NCO.h
 *
 *  Created on: Aug 17, 2019
 *      Author: bertold
 */

#ifndef NCOROTATE_H_
#define NCOROTATE_H_

#include <cmath>
#include <cstdint>

namespace NCO {

template <typename outType = float, typename phaseAccType = int32_t, int gain = 0> class NCORotate {
public:
	NCORotate(unsigned long sampleRate = 1):
			sampleRate_(sampleRate){

		v[0] = scaleFactor-1;
		v[1] = 0;

		r_c = scaleFactor;
		r_s = 0;
	}

	inline void fillBuffer(outType* sine, outType* cosine, unsigned int numSamples){
		for(unsigned int i=0; i<numSamples; i++){
			if(sine){
				*sine = v[1] >> gain;
				sine++;
			}

			if(cosine){
				*cosine = v[0] >> gain;
				cosine++;
			}

			phaseAccType amp = v[0]*v[0]+v[1]*v[1];

			auto tmp = (r_c*v[0] - r_s*v[1]);
			v[1]     = (r_s*v[0] + r_c*v[1]);
			v[0]=tmp;

			/* +1 to avoid overflow */
			phaseAccType s = 1 + (63*scaleFactor + amp/scaleFactor) / 64;

			v[0] /= s;
			v[1] /= s;
		}
	}

	inline void setFrequency(float freq){
		r_c = std::cos(freq/(float)sampleRate_ * 2.0 * M_PI) * (scaleFactor);
		r_s = std::sin(freq/(float)sampleRate_ * 2.0 * M_PI) * (scaleFactor);
	}

	float getAproxFrequency(){
		return std::atan2((scaleFactor-1)*r_s/scaleFactor, (scaleFactor-1)*r_c/scaleFactor) / 2.0 / M_PI * (float)sampleRate_;
	}

private:
	static const phaseAccType scaleFactor = (1<<(sizeof(phaseAccType)*4-1));

	unsigned int sampleRate_;

	phaseAccType v[2];
	phaseAccType r_c;
	phaseAccType r_s;
};

}

#endif /* NCOROTATE_H_ */

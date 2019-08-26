#include <cmath>
#include <type_traits>

#ifndef FMDISCRIMINATOR_H_
#define FMDISCRIMINATOR_H_

namespace Demod {

namespace FM {


enum FMDiscriminatorType{
	FM_DISC_ATAN,
	FM_DISC_DERIV,
};

template <FMDiscriminatorType type = FM_DISC_DERIV, typename InType = float, typename WorkType = float> class FMDiscriminator {
public:
	FMDiscriminator(){
	}

	inline void processBuffer(InType* I, InType* Q, InType* out, unsigned int numSamples){
		for(unsigned int k=0; k<numSamples; k++){
			WorkType i = I[k];
			WorkType q = Q[k];

			if(type == FM_DISC_ATAN){
				float at = std::atan2(i, q) / M_PI;
				float ato = at-oldAt;
				if(ato < -1){
					ato += 2;
				}
				*out = ato * (float)(1<<(sizeof(InType)*8));
				oldAt = at;

			}

			if(type == FM_DISC_DERIV){
				auto di = i -  delayI_[0];
				auto dq = q -  delayQ_[0];

				auto d = delayQ_[1] * di - delayI_[1] * dq;
				auto v = delayI_[1]*delayI_[1] + delayQ_[1]*delayQ_[1];

				if (!std::is_floating_point<InType>::value) {
					v /= 1<<(sizeof(InType)*8-2);
				}

				if(v){
					d /= v;
				}else{
					d = 0;
				}

				*out = d;

				delayI_[0] = delayI_[1];
				delayQ_[0] = delayQ_[1];
				delayI_[1] = i;
				delayQ_[1] = q;
			}

			out++;
		}
	}

private:
	WorkType delayI_[2]={};
	WorkType delayQ_[2]={};
	float oldAt = 0;
};

}

}


#endif /* FMDISCRIMINATOR_H_ */

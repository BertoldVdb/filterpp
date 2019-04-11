#ifndef IIRFILTER_H_
#define IIRFILTER_H_

#include <vector>

namespace Filter
{
namespace IIR
{

struct IIRCoefficients {
    std::vector<double> b;
    std::vector<double> a;
    double finalGain;
};

template <typename tapType = float, typename dataType = float, typename workType = float> class Filter
{
public:
    Filter(const IIRCoefficients& coeff):
        numCoeff_(coeff.a.size()),
        finalGain(coeff.finalGain)
    {
        a_ = vectorExtend<double, tapType>(coeff.a, 2, 1);
        b_ = vectorExtend<double, tapType>(coeff.b, 2, 1);
        buffer_ = std::vector<workType>(numCoeff_-1);
    }

    void filter(dataType* samplesIn, unsigned int numSamples, dataType* samplesOut)
    {
        for(unsigned int i=0; i<numSamples; i++) {
            workType feedback = samplesIn[i];
            workType out = 0;

            int l = bufferIndex_ + numCoeff_ - 1;
            for(unsigned int j=0; j<numCoeff_-1; j++) {
                feedback -= buffer_[j] * a_[l];
                out += buffer_[j] * b_[l];
                l--;
            }

            feedback /= a_[0];

            buffer_[bufferIndex_++] = feedback;
            if(bufferIndex_ >= numCoeff_ - 1 ) {
                bufferIndex_ = 0;
            }

            samplesOut[i] = (out + b_[0] * feedback) * finalGain;
        }
    }

private:
    unsigned int bufferIndex_ = 0;
    unsigned int numCoeff_;
    dataType finalGain;
    std::vector<tapType> a_ ,b_;
    std::vector<workType> buffer_;

    template <typename inType, typename outType> std::vector<outType> vectorExtend(std::vector<inType> in, unsigned int times = 2, unsigned int baseInc = 0, unsigned int endDec = 0)
    {
        std::vector<outType> out;

        out.insert(out.end(), in.begin(), in.end());

        for(unsigned int i=1; i<times; i++) {
            out.insert(out.end(), in.begin() + baseInc, in.end() - endDec);
        }

        return out;
    }
};

}
}



#endif /* IIRFILTER_H_ */

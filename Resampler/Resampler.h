#ifndef FILTER_RESAMPLER_RESAMPLER_H_
#define FILTER_RESAMPLER_RESAMPLER_H_

#include <memory>
#include <vector>

namespace Filter
{
namespace Resampler
{

template <typename T> struct ResamplerConfig {
    bool rational = true;
    double inputRate;
    double outputRate;

    unsigned int maxInterpolation = 32;
    unsigned int cicOrder = 6;
    unsigned int cicOversampleFactor = 4;

    double desiredRipple = -60;
    double transitionWidth = 1000;
    double cutoffFrequency = 5000;

    long cicBoost = 32768;

    bool filterInterpolate = false; /* This feature was promised, but is a serious computational burden, and normally not needed. */

    std::shared_ptr<std::vector<T>> taps;
};

}
}


#endif /* FILTER_RESAMPLER_RESAMPLER_H_ */

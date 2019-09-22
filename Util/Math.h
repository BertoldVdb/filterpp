#ifndef FILTERPP_UTIL_MATH_H_
#define FILTERPP_UTIL_MATH_H_

namespace Util {

/* This will be in std::numeric in C++17 :) */
template <typename T> T gcd(T v1, T v2)
{
    while (v2) {
        int tmp = v2;
        v2 = v1 % v2;
        v1 = tmp;
    }
    return v1;
}

}


#endif /* FILTERPP_UTIL_MATH_H_ */

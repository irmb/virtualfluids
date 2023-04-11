#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

#include "_deps/googletest-src/googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include <gmock/gmock.h>

inline auto RealEq = [](auto value) {
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value);
#else
    return testing::FloatEq(value);
#endif
};

inline auto RealNear = [](auto value, auto max_abs_error) {
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleNear(value, max_abs_error);
#else
    return testing::FloatNear(value, max_abs_error);
#endif
};

#endif

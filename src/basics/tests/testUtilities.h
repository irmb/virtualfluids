#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

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

namespace testingVF
{

__inline__ void captureStdOut()
{
    testing::internal::CaptureStdout();
}

__inline__ bool stdoutContainsWarning()
{
    std::string output = testing::internal::GetCapturedStdout();
    return output.find("warning") != std::string::npos;
}

} // namespace testingVF

#endif

#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

#include <gmock/gmock.h>
#include <basics/DataTypes.h>

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

namespace
{
__inline__ bool isEqualWithAccuracy(real number, real expected, real accuracy)
{
    return number > (expected - accuracy) && number < (expected + accuracy);
}
} // namespace

namespace testingVF
{

MATCHER_P2(RealNearForContainer, expectedContainer, accuracy, "")
{
    if (arg.size() != expectedContainer.size()) {
        std::cout << "The checked container does not have the same size as the expected container.\n" << std::endl;
        return false;
    }

    for (int i = 0; i < arg.size(); i++) {
        if (!isEqualWithAccuracy(arg[i], expectedContainer[i], accuracy)) {
            std::cout << "First mismatching element at index " << i << ": The actual element " << std::to_string(arg[i])
                      << " is not near the expected element " << std::to_string(expectedContainer[i])
                      << " (difference = " << std::to_string(expectedContainer[i] - arg[i]) << ").\n";
            return false;
        }
    }
    return true;
}

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

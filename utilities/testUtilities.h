#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

auto RealEq = [](auto value) {
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value);
#else
    return testing::FloatEq(value);
#endif
};

#endif
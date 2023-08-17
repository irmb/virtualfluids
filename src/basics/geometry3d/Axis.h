#ifndef AXIS
#define AXIS

#include <array>
#include <map>

enum Axis {
    x = 0,
    y = 1,
    z = 2,
};

const std::map<Axis, std::array<double, 3>> unitVectors{ { x, { 1, 0, 0 } },
                                                         { y, { 0, 1, 0 } },
                                                         { z, { 0, 0, 1 } } };

#endif
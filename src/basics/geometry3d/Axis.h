#ifndef AXIS
#define AXIS

#include <array>
#include <map>
#include <string>

enum Axis {
    x = 0,
    y = 1,
    z = 2,
};

namespace axis
{

const std::map<Axis, std::array<double, 3>> unitVectors{ { x, { 1, 0, 0 } },
                                                         { y, { 0, 1, 0 } },
                                                         { z, { 0, 0, 1 } } };

std::string to_string(Axis axis);

const std::array<Axis, 3> allAxes = { Axis::x, Axis::y, Axis::z };

} // namespace axis

#endif

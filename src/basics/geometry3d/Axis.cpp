#include "Axis.h"

std::string axis::to_string(Axis axis)
{
    switch (axis) {
        case x:
            return "x";
        case y:
            return "y";
        case z:
            return "z";
    }
    return "Axis not found.";
}
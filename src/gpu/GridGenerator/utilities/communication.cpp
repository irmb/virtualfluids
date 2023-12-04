#include "communication.h"

using namespace CommunicationDirections;

bool CommunicationDirections::isNegative(CommunicationDirection direction)
{
    return direction == CommunicationDirection::MX || direction == CommunicationDirection::MY ||
           direction == CommunicationDirection::MZ;
}

bool CommunicationDirections::isPositive(CommunicationDirection direction)
{
    return direction == CommunicationDirection::PX || direction == CommunicationDirection::PY ||
           direction == CommunicationDirection::PZ;
}

CommunicationDirection CommunicationDirections::getNegativeDirectionAlongAxis(Axis axis)
{
    switch (axis) {
        case Axis::x:
            return MX;
            break;
        case Axis::y:
            return MY;
            break;
        case Axis::z:
            return MZ;
            break;
        default:
            throw std::runtime_error("Unknown coordinate direction" + axis::to_string(axis));
    }
}

CommunicationDirection CommunicationDirections::getPositiveDirectionAlongAxis(Axis axis)
{
    switch (axis) {
        case Axis::x:
            return PX;
            break;
        case Axis::y:
            return PY;
            break;
        case Axis::z:
            return PZ;
            break;
        default:
            throw std::runtime_error("Unknown coordinate direction" + axis::to_string(axis));
    }
}
#ifndef DATA_POINT_H
#define DATA_POINT_H

#include <memory>

class DataPoint
{
public:
    static std::shared_ptr<DataPoint> getNewInstance(double x, double y);
    double getX();
    double getY();

private:
    DataPoint(double x, double y);
    DataPoint();

    double x;
    double y;
};
#endif
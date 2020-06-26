#ifndef CELL_H
#define CELL_H

#include "global.h"

#include "utilities/cuda/cudaDefines.h"

struct Point
{
    HOSTDEVICE Point() : x(0.0), y(0.0), z(0.0) {}
    HOSTDEVICE Point(real x, real y, real z) : x(x), y(y), z(z) {}
    real x, y, z;
};

class Cell
{
public:
    typedef Point* iterator;
    typedef const Point* const_iterator;

    HOSTDEVICE Cell(real startX, real startY, real startZ, real delta)
    {
        points = new Point[size];
        points[0] = Point(startX, startY, startZ); // 0,0,0
        points[1] = Point(startX + delta, startY, startZ); // 1,0,0
        points[2] = Point(startX, startY + delta, startZ); // 0,1,0
        points[3] = Point(startX + delta, startY + delta, startZ); // 1,1,0

        points[4] = Point(startX, startY, startZ + delta); // 0,0,1
        points[5] = Point(startX + delta, startY, startZ + delta); // 1,0,1
        points[6] = Point(startX, startY + delta, startZ + delta); // 0,1,1
        points[7] = Point(startX + delta, startY + delta, startZ + delta); // 1,1,1
    }

    HOSTDEVICE ~Cell()
    {
        delete[] points;
    }

    HOSTDEVICE iterator begin() { return &points[0]; }
    HOSTDEVICE const_iterator begin() const { return &points[0]; }
    HOSTDEVICE iterator end() { return &points[size]; }
    HOSTDEVICE const_iterator end() const { return &points[size]; }

private:
    Point* points;
    uint size = 8;

};



#endif
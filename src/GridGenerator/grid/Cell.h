#ifndef CELL_H
#define CELL_H


#include "core/DataTypes.h"


struct Point
{
    Point() : x(0.0), y(0.0), z(0.0) {}
    Point(real x, real y, real z) : x(x), y(y), z(z) {}
    real x, y, z;
};

class Cell
{
public:
    typedef Point* iterator;
    typedef const Point* const_iterator;

    Cell(real startX, real startY, real startZ, real delta)
    {
        points = new Point[size];
        points[0] = Point(startX, startY, startZ);
        points[1] = Point(startX + delta, startY, startZ);
        points[2] = Point(startX, startY + delta, startZ);
        points[3] = Point(startX + delta, startY + delta, startZ);

        points[4] = Point(startX, startY, startZ + delta);
        points[5] = Point(startX + delta, startY, startZ + delta);
        points[6] = Point(startX, startY + delta, startZ + delta);
        points[7] = Point(startX + delta, startY + delta, startZ + delta);
    }

    ~Cell()
    {
        delete[] points;
    }

    iterator begin() { return &points[0]; }
    const_iterator begin() const { return &points[0]; }
    iterator end() { return &points[size]; }
    const_iterator end() const { return &points[size]; }

private:
    Point* points;
    uint size = 8;

};








#endif

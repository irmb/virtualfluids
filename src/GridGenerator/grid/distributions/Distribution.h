#ifndef Distribution_H
#define Distribution_H

#include "GridGenerator/global.h"

#include <vector>
#include <string>

#define DIR_END_MAX 27


struct Direction
{
    HOSTDEVICE Direction()
    {
        dir[0] = 0;
        dir[1] = 0;
        dir[2] = 0;
    }

    HOSTDEVICE Direction(int dx, int dy, int dz)
    {
        dir[0] = dx;
        dir[1] = dy;
        dir[2] = dz;
    }

    HOSTDEVICE int operator[](uint dir) const
    {
        if (dir < 3)
            return this->dir[dir];
        return -99;
    }
private:
    int dir[3];
};

struct Distribution
{
    typedef Direction* iterator;
    typedef const Direction* const_iterator;

    real* f;
    int *dirs;
    Direction* directions;
    int dir_start;
    int dir_end;
    const char* name;
    unsigned long fSize;

    void setSize(uint size) {
        fSize = size * (dir_end + 1);
    }

    HOSTDEVICE iterator begin() { return &directions[0]; }
    HOSTDEVICE const_iterator begin() const { return &directions[0]; }
    HOSTDEVICE iterator end() { return &directions[dir_end + 1]; }
    HOSTDEVICE const_iterator end() const { return &directions[dir_end + 1]; }
};

class Grid;

class VF_PUBLIC DistributionHelper
{
public:
    static Distribution getDistribution7();
    static Distribution getDistribution13();
    static Distribution getDistribution19();
    static Distribution getDistribution27();

    static Distribution getDistribution(std::string name);

public:
    static std::vector<std::vector<real> > getQsWithoutRowsWithOnlyZeroValues(const Grid &grid, const Distribution &d);
    static std::vector<std::vector<real> > getAllQsOnFluidNodes(const Grid &grid, const Distribution &d);
    static int getNeighborNodeIndexInGivenDirection(const Distribution &d, const Grid &grid, const int node, const int dir_index);
    static std::vector<std::vector<real> > getVectorWithoutRowsWithOnlyZeroValues(std::vector<std::vector<real> > qs);
    static void printQs(std::vector<std::vector<real> > qs, int decimalPlaces);
};

#endif

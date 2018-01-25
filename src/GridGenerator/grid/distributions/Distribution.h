#ifndef Distribution_H
#define Distribution_H

#define DIR_END_MAX 27

#include "GridGenerator/global.h"


#include <vector>
#include <string>


struct Distribution
{
    real* f;
    int *dirs;
    int dir_start;
    int dir_end;
    const char* name;
    unsigned long fSize;

    void setSize(uint size) {
        fSize = size * (dir_end + 1);
    }
};

struct Grid;

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

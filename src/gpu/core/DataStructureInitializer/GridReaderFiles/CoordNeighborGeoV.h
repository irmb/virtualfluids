#ifndef CoordNeighborGeoV_H
#define CoordNeighborGeoV_H

#include <vector>
#include <string>
#include <fstream>
#include "Calculation/Calculation.h"


class CoordNeighborGeoV 
{
protected:
    std::ifstream file; 
    unsigned int maxLevel; 
    std::vector<unsigned int> levelSizes;
    std::vector< std::vector<unsigned int> > neighbors;
    std::vector< std::vector< real> > coordinates; 

public:
    CoordNeighborGeoV();
    CoordNeighborGeoV(std::string ad, bool binaer, bool coord);
    ~CoordNeighborGeoV(void);

    void init(bool coord);
    void init_Binary(bool coord);

    unsigned int getLevel();
    unsigned int getSize(unsigned int level);
    std::vector<unsigned int>getVec(unsigned int level);
    void setVec(unsigned int level, std::vector<unsigned int> vec);

    void initalNeighbors(unsigned int *int_ptr, unsigned int level ) const;
    void initalCoords(real *int_ptr, unsigned int level ) const;

protected:
    void skipSpace();
    void readLevelSize(unsigned int level);
    void readNeighbors(unsigned int level);
    void readCoordinates(unsigned int level);
    void resizeVectors();
    void readLevel();

};

#endif

#ifndef OFFSET_SCALE_H
#define OFFSET_SCALE_H

#include <vector>
#include <string>

#include "CoordNeighborGeoV.h"

class OffsetScale :
	public CoordNeighborGeoV
{
private:
	std::vector<real>vec1D_dataOffset; //das Feld mit Werten, tempor�r zum F�llen von vec2D
    std::vector<std::vector<real> >vec2D_dataOffset; //alle Felder mit Werten gegliedert nach Level
    std::vector<std::vector<unsigned int> > scale;
public:
	OffsetScale(std::string ad, bool off);
	~OffsetScale(void);
	void init();

    void initScale(unsigned int* data, unsigned int level);

    void initOffset();
	void initArrayOffset(real *x_ptr,real *y_ptr,real *z_ptr, unsigned int level);
};

#endif

#pragma once
#include "CoordNeighborGeoV.h"
class OffsetScale :
	public CoordNeighborGeoV
{
private:
	vector<doubflo>vec1D_dataOffset; //das Feld mit Werten, temporär zum Füllen von vec2D
	vector<vector<doubflo> >vec2D_dataOffset; //alle Felder mit Werten gegliedert nach Level
public:
	OffsetScale(string ad, bool off);
	~OffsetScale(void);
	void init();

	void initOffset();
	void initArrayOffset(doubflo *x_ptr,doubflo *y_ptr,doubflo *z_ptr, unsigned int level);
};


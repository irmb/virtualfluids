#ifndef GRID_INFORMATION_STRUCT_H
#define GRID_INFORMATION_STRUCT_H

#include <string>

struct GridInformationStruct
{
	unsigned int numberOfGridLevels;
	unsigned int maxLevel;
	std::string gridPath;
	double lx;
	double lz;
};
#endif
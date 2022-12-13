#ifndef DISTRIBUTION_INSPECTOR_H
#define DISTRIBUTION_INSPECTOR_H

#include "Parameter/Parameter.h"


class DistributionDebugInspector
{
public:
	DistributionDebugInspector(uint _inspectionLevel, real _minX, real _maxX, real _minY, real _maxY, real _minZ, real _maxZ, std::string _tag):
    inspectionLevel(_inspectionLevel),
    minX(_minX),
    maxX(_maxX),
    minY(_minY),
    maxY(_maxY),
    minZ(_minZ),
    maxZ(_maxZ),
    tag(_tag)
    {};
	
    ~DistributionDebugInspector(){}

    void inspect(std::shared_ptr<Parameter> para, uint level, uint t);


private:
uint inspectionLevel;
real minX;
real maxX;
real minY;
real maxY;
real minZ;
real maxZ;
std::string tag;

};

#endif
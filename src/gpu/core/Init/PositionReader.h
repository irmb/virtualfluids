#ifndef POSITION_READER_H
#define POSITION_READER_H

class Parameter;

class PositionReader
{
public:
   static void readFilePropellerCylinderForAlloc(Parameter* para);
   static void readFilePropellerCylinder(Parameter* para);
   static void definePropellerQs(Parameter* para);
   static void readMeasurePoints(Parameter* para);
};

#endif
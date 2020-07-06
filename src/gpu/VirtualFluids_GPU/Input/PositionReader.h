#ifndef POSITION_READER_H
#define POSITION_READER_H

#include <string>
#include "Parameter/Parameter.h"

class PositionReader
{
public:
   PositionReader(){}
   ~PositionReader(){}
   static void readFileForAlloc(const std::string fileName, Parameter* para);
   static void readFile(const std::string fileName, std::string Type, Parameter* para);
   static void readFileInterfaceForAlloc(const std::string fileName, std::string Type, Parameter* para);
   static void readFileInterface(const std::string fileName, std::string Type, Parameter* para);
   static void readFileInterfaceOffsetForAlloc(const std::string fileName, std::string Type, Parameter* para);
   static void readFileInterfaceOffset(const std::string fileName, std::string Type, Parameter* para);
   static void readFileNoSlipBcForAlloc(const std::string fileName, Parameter* para);
   static void readFileNoSlipBcQreadForAlloc(const std::string fileName, Parameter* para);
   static void readFileNoSlipBcPos(const std::string fileName, Parameter* para);
   static void readFileNoSlipBcQs(const std::string fileName, Parameter* para);
   static void readFileNoSlipBcValue(const std::string fileName, Parameter* para);
   static void findQs(Parameter* para);
   static void readFileSlipBcForAlloc(const std::string fileName, Parameter* para);
   static void readFileSlipBcQreadForAlloc(const std::string fileName, Parameter* para);
   static void readFileSlipBcPos(const std::string fileName, Parameter* para);
   static void readFileSlipBcQs(const std::string fileName, Parameter* para);
   static void readFileSlipBcValue(const std::string fileName, Parameter* para);
   static void findSlipQs(Parameter* para);
   static void readFilePressBcForAlloc(const std::string fileName, Parameter* para);
   static void readFilePressBcQreadForAlloc(const std::string fileName, Parameter* para);
   static void readFilePressBcPos(const std::string fileName, Parameter* para);
   static void readFilePressBcQs(const std::string fileName, Parameter* para);
   static void readFilePressBcValue(const std::string fileName, Parameter* para);
   static void findPressQs(Parameter* para);
   static void readFilePropellerCylinderForAlloc(Parameter* para);
   static void readFilePropellerCylinder(Parameter* para);
   static void definePropellerQs(Parameter* para);
   static void readMeasurePoints(Parameter* para);
protected:
private:
};

#endif
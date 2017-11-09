#ifndef KFULL_READER_H
#define KFULL_READER_H

#include <string>
#include "Parameter/Parameter.h"

class kFullReader
{
public:
   kFullReader(){}
   ~kFullReader(){}
   static void readFile(const std::string fileName, Parameter* para);
   static void readFileForAlloc(const std::string fileName, Parameter* para);
   static void readGeoFull(const std::string fileName, Parameter* para);
protected:
private:
};

#endif
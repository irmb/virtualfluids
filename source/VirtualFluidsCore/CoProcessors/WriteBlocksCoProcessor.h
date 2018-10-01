/*
*  WriteBlocksCoProcessor.h
*
*  Created on: 24.09.2012
*  Author: K. Kucher
*/

#ifndef WriteBlocksCoProcessor_H_
#define WriteBlocksCoProcessor_H_

#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;

class WriteBlocksCoProcessor: public CoProcessor 
{
public:
   WriteBlocksCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<Communicator> comm);
   virtual ~WriteBlocksCoProcessor();

   void process(double step) override;

protected:
   void collectData(double step);

   std::string path;
   WbWriter* writer;
   SPtr<Communicator>  comm;
};


#endif 

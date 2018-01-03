/*
*  BlocksCoProcessor.h
*
*  Created on: 24.09.2012
*  Author: K. Kucher
*/

#ifndef BlocksCoProcessor_H_
#define BlocksCoProcessor_H_

#include <memory>
#include <string>

#include "CoProcessor.h"

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;

class WriteBlocksCoProcessor;
typedef std::shared_ptr<WriteBlocksCoProcessor> WriteBlocksCoProcessorPtr;

class WriteBlocksCoProcessor: public CoProcessor 
{
public:
   WriteBlocksCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, const std::string& path, WbWriter* const writer, std::shared_ptr<Communicator> comm);
   virtual ~WriteBlocksCoProcessor();

   void process(double step) override;

protected:
   void collectData(double step);

   std::string path;
   WbWriter* writer;
   std::shared_ptr<Communicator>  comm;
};


#endif 

/*
 *  EmergencyExitCoProcessor.h
 *
 *  Created on: 05.10.2012
 *  Author: K. Kucher
 */

#ifndef EmergencyExitCoProcessor_H
#define EmergencyExitCoProcessor_H

#include "CoProcessor.h"
#include "Communicator.h"
#include "RestartCoProcessor.h"

#include <boost/shared_ptr.hpp>
class EmergencyExitCoProcessor;
typedef boost::shared_ptr<EmergencyExitCoProcessor> EmergencyExitCoProcessorPtr;

class EmergencyExitCoProcessor: public CoProcessor 
{
public:
	EmergencyExitCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
                              const std::string& path, RestartCoProcessorPtr rp,
                              CommunicatorPtr comm);
	virtual ~EmergencyExitCoProcessor();
	void process(double step);
protected:
	void collectData(double step);
   void writeMetafile(int status);
   bool readMetafile();
   void checkMetafile();
private:
   std::string path;
   CommunicatorPtr comm;
   RestartCoProcessorPtr rp;
   std::string metafile;
};


#endif /* EmergencyExitCoProcessor_H */

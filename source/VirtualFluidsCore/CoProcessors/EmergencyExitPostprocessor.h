/*
 *  EmergencyExitPostprocessor.h
 *
 *  Created on: 05.10.2012
 *  Author: K. Kucher
 */

#ifndef EmergencyExitPostprocessor_H
#define EmergencyExitPostprocessor_H

#include "Postprocessor.h"
#include "Communicator.h"
#include "RestartPostprocessor.h"

#include <boost/shared_ptr.hpp>
class EmergencyExitPostprocessor;
typedef boost::shared_ptr<EmergencyExitPostprocessor> EmergencyExitPostprocessorPtr;

class EmergencyExitPostprocessor: public Postprocessor 
{
public:
	EmergencyExitPostprocessor(Grid3DPtr grid, UbSchedulerPtr s,
                              const std::string& path, RestartPostprocessorPtr rp,
                              CommunicatorPtr comm);
	virtual ~EmergencyExitPostprocessor();
	void update(double step);
protected:
	void collectPostprocessData(double step);
   void writeMetafile(int status);
   bool readMetafile();
   void checkMetafile();
private:
   std::string path;
   CommunicatorPtr comm;
   RestartPostprocessorPtr rp;
   std::string metafile;
};


#endif /* EmergencyExitPostprocessor_H */

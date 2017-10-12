//! \file QCriterionCoProcessor.h
//!  \brief Created on: 25.08.2013
//!  \author: Sonja Uphoff


#ifndef QCriterionCoProcessor_H
#define QCriterionCoProcessor_H

#include "CoProcessor.h"
#include "Grid3D.h"
#include "Block3D.h"

#include "Communicator.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class QCriterionCoProcessor;
typedef boost::shared_ptr<QCriterionCoProcessor> QCriterionCoProcessorPtr;

//! \brief  Computes the value Q with which vortices can be visualized as isocontours to Q=0, writes to .vtk, For uniform, serial setups only!
//! \details writes at given time intervals specified in scheduler (s)  
//!          Processing: paraview, take isolines of entry for Q-criterion vortex detection 
//!			 Q-Criterion: Visualize Vorteces as regions where Vorticity is larger than strain rate (Hunt, 1988)
//! \author  Sonja Uphoff 

class QCriterionCoProcessor : public CoProcessor
{
public:
	QCriterionCoProcessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer, 
		UbSchedulerPtr s, CommunicatorPtr comm);
	//! Make update if timestep is write-timestep specified in UbSchedulerPtr s
	void process(double step); 
protected:
	//! Prepare data and write in .vtk file
	void collectData(double step);
	//! Q is computed for all points in a block. Data for writing is added to data and cell vectors. 
	void addData(const Block3DPtr block);
	//! After writing to .vtk-file, all vectors are reset 
	void clearData();
	//! Computes macroscopic velocities 
	void computeVelocity(LBMReal* f, LBMReal* v);
	//! Computes average and RMS values of macroscopic quantities 
	void getNeighborVelocities(int offx, int offy, int offz, int ix1, int ix2, int ix3,const Block3DPtr block, LBMReal* vE,LBMReal* vW);

private:
	void init();
	std::vector<UbTupleFloat3> nodes;
	std::vector<UbTupleInt8> cells;
	std::vector<std::string> datanames; //only one entry for QKrit-CoProcessor: Q
	std::vector<std::vector<double> > data; 
	std::vector<std::vector<Block3DPtr> > blockVector;
	int minInitLevel; //go through all levels for block vector of current process from minInitLevel to maxInitLevel
	int maxInitLevel;
	int gridRank;     //comm-Rank des aktuellen prozesses
	std::string path;
	WbWriter* writer;
	CommunicatorPtr comm;
	enum Values{xdir = 0, ydir = 1, zdir = 2};  	//labels for the different components
};
#endif



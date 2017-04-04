#ifndef D3Q27ShearStressCoProcessor_H
#define D3Q27ShearStressCoProcessor_H

#include "CoProcessor.h"
#include "Communicator.h"
#include "D3Q27Interactor.h"
#include "InterpolationProcessor.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class ShearStressCoProcessor;
typedef boost::shared_ptr<ShearStressCoProcessor> ShearStressCoProcessorPtr;

//! \brief  Computes the shear stress and y plus values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s) and resets according to scheduler (rs).  
//!          Take root to obtain  during post processing (paraview).   
//! \author  K. Kucher, S. Uphoff, M. Geier, E. Goraki Fard  

class ShearStressCoProcessor: public CoProcessor 
{
public:
   //! Defoult constructor
   ShearStressCoProcessor(){}
   //! Constructor
   ShearStressCoProcessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer, 
      UbSchedulerPtr s, UbSchedulerPtr rs);
   virtual ~ShearStressCoProcessor();             
   void process(double step); 
   void addInteractor(D3Q27InteractorPtr interactor);
protected:
   //! Computes average and shear stress values of macroscopic quantities 
   void calculateShearStress(double timeStep);
   //! Prepare data and write in .vtk file
   void collectData(double step);
   //! Reset data
   void resetData(double step);
   //! prepare data
   void addData();
   void clearData();
   void reset(double step);
   void findPlane(int ix1,int ix2,int ix3,Grid3DPtr grid,Block3DPtr block,double &A,double &B,double &C,double &D,double &ii);
   bool checkUndefindedNodes( BCArray3D& bcArray,int ix1,int ix2,int ix3);
   void initDistance();
private:
   std::vector<UbTupleFloat3> nodes;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;
   std::string path;
   std::vector<D3Q27InteractorPtr> interactors;
   std::vector<double> normals;
   int gridRank;
   WbWriter* writer;
   UbSchedulerPtr Resetscheduler;  //additional scheduler to restart averaging after a given interval
   int minInitLevel; //min init level
   int maxInitLevel;
   std::vector<std::vector<Block3DPtr> > blockVector;
   enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvSxx = 3, AvSyy = 4, AvSzz = 5, AvSxy = 6, AvSyz = 7, AvSxz = 8, normalX1 = 9, normalX2 = 10, normalX3 = 11, normalq = 12,numberOfPoint=13}; 

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<CoProcessor>(*this);
      ar & path;
      ar & normals;
      ar & interactors;
      ar & blockVector;
      ar & minInitLevel;
      ar & maxInitLevel;
      ar & gridRank;
      ar & writer;
   }
};


#endif /* D3Q27ShearStressCoProcessor_H */

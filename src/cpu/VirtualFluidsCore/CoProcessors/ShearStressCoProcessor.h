#ifndef D3Q27ShearStressCoProcessor_H
#define D3Q27ShearStressCoProcessor_H

#include <PointerDefinitions.h>
#include <vector>
#include <string>

#include <basics/utilities/UbTuple.h>

#include "CoProcessor.h"

class Block3D;
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class BCArray3D;
class WbWriter;

//! \brief  Computes the shear stress and y plus values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s) and resets according to scheduler (rs).  
//!          Take root to obtain  during post processing (paraview).   
//! \author  K. Kucher, S. Uphoff, M. Geier, E. Goraki Fard  

class ShearStressCoProcessor: public CoProcessor 
{
public:
   //! Default constructor
   ShearStressCoProcessor()= default;
   //! Constructor
   ShearStressCoProcessor(SPtr<Grid3D> grid, const std::string& path, WbWriter* const writer,
       SPtr<UbScheduler> s, SPtr<UbScheduler> rs);
   ~ShearStressCoProcessor() override; 
    
   void process(double step) override; 

   void addInteractor(SPtr<D3Q27Interactor> interactor);
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
   void findPlane(int ix1,int ix2,int ix3, SPtr<Grid3D> grid, SPtr<Block3D> block,double &A,double &B,double &C,double &D,double &ii);
   bool checkUndefindedNodes(SPtr<BCArray3D> bcArray,int ix1,int ix2,int ix3);
   void initDistance();

private:
   std::vector<UbTupleFloat3> nodes;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;
   std::string path;
   std::vector<SPtr<D3Q27Interactor> > interactors;
   std::vector<double> normals;
   int gridRank;
   WbWriter* writer;
   SPtr<UbScheduler> Resetscheduler;  //additional scheduler to restart averaging after a given interval
   int minInitLevel; //min init level
   int maxInitLevel;
   std::vector<std::vector<SPtr<Block3D> > > blockVector;
   enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvSxx = 3, AvSyy = 4, AvSzz = 5, AvSxy = 6, AvSyz = 7, AvSxz = 8, normalX1 = 9, normalX2 = 10, normalX3 = 11, normalq = 12,numberOfPoint=13}; 
};


#endif /* D3Q27ShearStressCoProcessor_H */

#ifdef VF_VTK

#ifndef InSituVTKCoProcessor_h__
#define InSituVTKCoProcessor_h__

#include <CoProcessor.h>
#include <Grid3D.h>
#include <LBMUnitConverter.h>

#include <string>

//VTK headers
#include <vtkSocketCommunicator.h>
#include <vtkSocketController.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>

class InSituVTKCoProcessor : public CoProcessor
{
public:
   InSituVTKCoProcessor();
   InSituVTKCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& configFile, SPtr<LBMUnitConverter> conv);
   virtual ~InSituVTKCoProcessor(); 
   void process(double step);
protected:
   void collectData(double step);
   void addData(SPtr<Block3D> block);
   void readConfigFile(const std::string& configFile);

   //void clearData();
private:
   std::string path;
   SPtr<LBMUnitConverter> conv;
   std::vector<std::vector<SPtr<Block3D>> > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;
   vtkSmartPointer<vtkSocketCommunicator> comm;
   vtkSmartPointer<vtkSocketController>   contr;
   vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
   vtkSmartPointer<vtkPoints> points;
   vtkSmartPointer<vtkDoubleArray> arrays[5];
   int wPort;
   std::string wHostname;
   std::string wIP;
};

#endif // InSituVTKCoProcessor_h__

#endif


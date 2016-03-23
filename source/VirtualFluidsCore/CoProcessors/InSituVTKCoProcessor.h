#ifdef VF_VTK

#ifndef InSituVTKPostprocessor_h__
#define InSituVTKPostprocessor_h__

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
   InSituVTKCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& configFile, LBMUnitConverterPtr conv);
   virtual ~InSituVTKCoProcessor(); 
   void process(double step);
protected:
   void collectData(double step);
   void addData(Block3DPtr block);
   void readConfigFile(const std::string& configFile);

   //void clearData();
private:
   std::string path;
   LBMUnitConverterPtr conv;
   std::vector<std::vector<Block3DPtr> > blockVector;
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

#endif // InSituVTKPostprocessor_h__

#endif


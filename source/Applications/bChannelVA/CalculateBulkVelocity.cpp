#include "Postprocessing.h"
#include <vtkPlane.h>
#include <vtkClipDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>

void CalculateBulkVelocity(std::string dataName, double origin[3])
{
   vtkSmartPointer<vtkDataSet> dataSet(ReadDataSet(dataName.c_str()));

   //vtkPlane 
   vtkSmartPointer<vtkPlane> plane = vtkPlane::New();
   plane->SetOrigin(origin[0], origin[1], origin[2]);
   plane->SetNormal(0.0, 0.0, 1.0);

   //Clip
   vtkSmartPointer<vtkClipDataSet> clipper = vtkClipDataSet::New();
   clipper->SetInputData(dataSet);
   clipper->SetClipFunction(plane);
   clipper->SetValue(0.5);
   clipper->Update();

   double numPoints = clipper->GetOutput()->GetNumberOfPoints();
   vtkSmartPointer<vtkDataArray> vxArray = clipper->GetOutput()->GetPointData()->GetArray("vx");
   double vxAvSumm = 0;
   double vxAv = 0;
   std::cout << "numPints = "<<numPoints<<std::endl;

   for (int i = 0; i < (int)numPoints; i++)
   {
      vxAvSumm += vxArray->GetTuple1(i);
      //std::cout << "Point = "<<i<<std::endl;
   }

   if(numPoints > 0)
   {
      vxAv = vxAvSumm/numPoints;
   }


   std::cout << "vxAvSumm = "<<vxAvSumm<<std::endl;
   std::cout << "vxAv = "<<vxAv<<std::endl;
}


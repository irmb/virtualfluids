#include "ReadDataSet.h"
#include <vtkXMLReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLRectilinearGridReader.h>
//#include <vtkXMLHyperOctreeReader.h>
#include <vtkXMLCompositeDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
//#include <vtkHyperOctree.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtksys/SystemTools.hxx>

vtkDataSet* ReadDataSet(std::string fileName)
{
   std::string extension =
      vtksys::SystemTools::GetFilenameLastExtension(fileName);
   // Dispatch based on the file extension
   if (extension == ".vtu")
   {
      return ReadAnXMLFile<vtkXMLUnstructuredGridReader> (fileName);
   }
   else if (extension == ".vtp")
   {
      return ReadAnXMLFile<vtkXMLPolyDataReader> (fileName);
   }
   else if (extension == ".vts")
   {
      return ReadAnXMLFile<vtkXMLStructuredGridReader> (fileName);
   }
   else if (extension == ".vtr")
   {
      return ReadAnXMLFile<vtkXMLRectilinearGridReader> (fileName);
   }
   else if (extension == ".vti")
   {
      return ReadAnXMLFile<vtkXMLImageDataReader> (fileName);
   }
   //else if (extension == ".vto")
   //{
   //   return ReadAnXMLFile<vtkXMLHyperOctreeReader> (fileName);
   //}
   else if (extension == ".vtk")
   {
      return ReadAnXMLFile<vtkDataSetReader> (fileName);
   }
   else if (extension == ".pvtu")
   {
      return ReadAnXMLFile<vtkXMLPUnstructuredGridReader> (fileName);
   }
   else if (extension == ".pvtp")
   {
      return ReadAnXMLFile<vtkXMLPPolyDataReader> (fileName);
   }
   else
   {
      return NULL;
   }
}


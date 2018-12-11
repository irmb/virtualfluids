#include <string>
#include <vtkDataSet.h>

vtkDataSet* ReadDataSet(std::string fileName);

//////////////////////////////////////////////////////////////////////////
template<class TReader> vtkDataSet* ReadAnXMLFile(std::string fileName)
{
   vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
   reader->SetFileName(fileName.c_str());
   reader->Update();
   reader->GetOutput()->Register(reader);
   return vtkDataSet::SafeDownCast(reader->GetOutput());
}
//////////////////////////////////////////////////////////////////////////
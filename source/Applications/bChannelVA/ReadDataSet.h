#include <string>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>

vtkDataSet* ReadDataSet(std::string fileName);

//////////////////////////////////////////////////////////////////////////
template<class TReader> 
vtkDataSet* ReadAnXMLFile(std::string fileName)
{
   TReader* reader = TReader::New();
   reader->SetFileName(fileName.c_str());
   reader->Update();
   reader->GetOutput()->Register(reader);
   return vtkDataSet::SafeDownCast(reader->GetOutput());
}
//////////////////////////////////////////////////////////////////////////
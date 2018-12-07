#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>

using namespace std;

void StlToVtu(string stlfile, string vtifile, double spacing[3])
{
   //read STL-file
   vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
   reader->SetFileName(stlfile.c_str());
   reader->Update();

   vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
   reader->Update();

   //create image data
   vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    
   double bounds[6];
   pd->GetBounds(bounds);
   whiteImage->SetSpacing(spacing);

   // compute dimensions
   int dim[3];
   for (int i = 0; i < 3; i++)
   {
      dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
   }
   whiteImage->SetDimensions(dim);
   whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

   double origin[3];
   origin[0] = bounds[0] + spacing[0] / 2;
   origin[1] = bounds[2] + spacing[1] / 2;
   origin[2] = bounds[4] + spacing[2] / 2;
   whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
   whiteImage->SetScalarTypeToUnsignedChar();
   whiteImage->AllocateScalars();
#else
   whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
   // fill the image with foreground voxels:
   unsigned char inval = 255;
   unsigned char outval = 0;
   vtkIdType count = whiteImage->GetNumberOfPoints();
   for (vtkIdType i = 0; i < count; ++i)
   {
      whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
   }

   // polygonal data --> image stencil:
   vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
   pol2stenc->SetInput(pd);
#else
   pol2stenc->SetInputData(pd);
#endif
   pol2stenc->SetOutputOrigin(origin);
   pol2stenc->SetOutputSpacing(spacing);
   pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
   pol2stenc->Update();

   // cut the corresponding white image and set the background:
   vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
   imgstenc->SetInput(whiteImage);
   imgstenc->SetStencil(pol2stenc->GetOutput());
#else
   imgstenc->SetInputData(whiteImage);
   imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
   imgstenc->ReverseStencilOff();
   imgstenc->SetBackgroundValue(outval);
   imgstenc->Update();

   vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
   writer->SetInputData(imgstenc->GetOutput());
   writer->SetFileName(vtifile.c_str());
   writer->SetDataModeToAscii();
   writer->Update();
}

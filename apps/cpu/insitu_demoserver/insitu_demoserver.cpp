#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include "vtkActor.h"
#include "vtkContourFilter.h"
#include "vtkDataSetMapper.h"
#include "vtkDebugLeaks.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRectilinearGrid.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderWindow.h"
#include "vtkOpenGLRenderer.h"
#include "vtkSocketCommunicator.h"
#include "vtkSocketController.h"
#include "vtkStructuredGrid.h"
#include "vtkImageData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCamera.h"
#include "vtkImageActor.h"
#include <vtkXMLUnstructuredGridWriter.h>
#include "vtkRenderWindowInteractor.h"
#include "vtkOpenGLActor.h"
#include "vtkSmartPointer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkPlane.h>
#include <vtkCutter.h>

#include <boost/thread.hpp>

#define VTK_CREATE(type, name) \
   vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


static const int scMsgLength = 10;

static void CleanUp(vtkSmartPointer<vtkSocketCommunicator> vtkNotUsed(comm),
                    vtkSmartPointer<vtkSocketController> vtkNotUsed(contr))
{
   // This will close the connection as well as delete
   // the communicator
   // Deleting no longer necessary with smart pointers.
   //   comm->Delete();
   //   contr->Delete();
}

using namespace std;

vtkSmartPointer<vtkSocketController> contr;
vtkSmartPointer<vtkSocketCommunicator> comm;

void receive(vtkSmartPointer<vtkUnstructuredGrid> ugrid, vtkSmartPointer<vtkDataSetMapper> umapper, vtkSmartPointer<vtkRenderWindow> renWin)
{
   int step;
   while (true)
   {
      if (!comm->Receive(&step, 1, 1, 11))
      {
         cerr << "Server error: Error receiving data." << endl;
         CleanUp(comm, contr);
         return;
      }

      cout << "step: "<<step<<"\n";

      if (!comm->Receive(ugrid, 1, 9))
      {
         cerr << "Client error: Error receiving data." << endl;
         CleanUp(comm, contr);
         return;
      }
      double range[2];
      ugrid->GetPointData()->GetArray("Vx")->GetRange(range);
      umapper->SetScalarRange(range);
      umapper->Update();
      //renWin->Render();
   }
}

////////////////////////////////////////////////////////////////////////
void server()
{
   try
   {
      contr = vtkSmartPointer<vtkSocketController>::New();
      contr->Initialize();

      comm = vtkSmartPointer<vtkSocketCommunicator>::New();

      string hostname = "localhost";
      int port=11111;

      // Establish connection
      if (!comm->WaitForConnection(port))
      {
         cerr << "Server error: Wait timed out or could not initialize socket." << endl;
         return;
      }

      // Test receiving vtkDataObject
      VTK_CREATE(vtkUnstructuredGrid, ugrid);

      int step;

      if (!comm->Receive(&step, 1, 1, 11))
      {
         cerr << "Server error: Error receiving data." << endl;
         CleanUp(comm, contr);
         return;
      }

      cout << "step: "<<step<<"\n";

      if (!comm->Receive(ugrid, 1, 9))
      {
         cerr << "Client error: Error receiving data." << endl;
         CleanUp(comm, contr);
         return;
      }

      vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      writer->SetInput(ugrid);
      writer->SetFileName("test.vtu");
      writer->SetDataModeToAscii();
      writer->Update();

      //vtkPlane
      vtkSmartPointer<vtkPlane> plane = vtkPlane::New();
      plane->SetNormal(0.0, 1.0, 0.0);
      plane->SetOrigin(40, 19.5, 19.5);

      //Cut
      vtkSmartPointer<vtkCutter> planeCut = vtkCutter::New();
      planeCut->SetInput(ugrid);
      planeCut->SetCutFunction(plane);
      planeCut->Update();

      VTK_CREATE(vtkDataSetMapper, umapper);
      //umapper->SetInput(planeCut->GetOutput());
      umapper->SetInput(ugrid);

      umapper->SetScalarModeToUsePointFieldData();
      umapper->SetColorModeToMapScalars();
      umapper->ScalarVisibilityOn();
      double range[2];
      //planeCut->GetOutput()->GetPointData()->GetArray("Vx")->GetRange(range);
      ugrid->GetPointData()->GetArray("Vx")->GetRange(range);
      umapper->SetScalarRange(range);
      umapper->SelectColorArray("Vx");

      VTK_CREATE(vtkActor, uactor);
      uactor->SetMapper(umapper);

      VTK_CREATE(vtkRenderer, ren);
      ren->AddActor(uactor);
      ren->SetBackground( 0.1, 0.2, 0.4 );

      VTK_CREATE(vtkRenderWindow, renWin);
      renWin->SetSize(1024,800);
      renWin->AddRenderer(ren);

      //while (true)
      //{
      //   if (!comm->Receive(&step, 1, 1, 11))
      //   {
      //      cerr << "Server error: Error receiving data." << endl;
      //      CleanUp(comm, contr);
      //      return;
      //   }

      //   cout << "step: "<<step<<"\n";

      //   if (!comm->Receive(ugrid, 1, 9))
      //   {
      //      cerr << "Client error: Error receiving data." << endl;
      //      CleanUp(comm, contr);
      //      return;
      //   }

      //   //writer->Update();
      //   
      //   planeCut->Update();
      //   planeCut->GetOutput()->GetPointData()->GetArray("Vx")->GetRange(range);
      //   umapper->SetScalarRange(range);
      //   umapper->Update();
      //   renWin->Render();
      //}

      boost::thread t(boost::bind( &receive, ugrid, umapper,renWin));

      vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
      iren->SetRenderWindow(renWin);

      vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
      iren->SetInteractorStyle(style);

      iren->Initialize();
      iren->Start();

      iren->Delete();
      style->Delete();

      CleanUp(comm, contr);

   }
   catch(std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch(std::string& s)
   {
      cerr << s << endl;
   }
   catch(...)
   {
      cerr << "unknown exception" << endl;
   }

}

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   server();
}


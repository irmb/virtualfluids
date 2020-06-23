#include <iostream>
#include <string>

#include "numerics/geometry3d/CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "GenBlocksGridVisitor.h"
#include "numerics/geometry3d/GbSystem3D.h"
#include "numerics/geometry3d/GbCuboid3D.h"
#include "numerics/geometry3d/GbCylinder3D.h"
#include <numerics/geometry3d/GbSphere3D.h>
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "RefineInterGbObjectsVisitor.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "SetKernelBlockVisitor.h"
#include "LBMKernelETD3Q27Cascaded.h"
#include "D3Q27MacroscopicQuantitiesPostprocessor.h"
#include "MPICommunicator.h"
#include "D3Q27ETBCProcessor.h"
#include "SimulationParameters.h"
#include "D3Q27SetUndefinedNodesBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include "D3Q27SetConnectorsBlockVisitor.h"
#include "NullCommunicator.h"
#include "D3Q27ETInitDistributionsBlockVisitor.h"
#include "CalculationManager.h"
#include "PQueuePartitioningGridVisitor.h"
#include "MetisPartitioningGridVisitor.h"
#include "D3Q27Interactor.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27VelocityBCAdapter.h"
#include "D3Q27DensityBCAdapter.h"
#include "D3Q27BoundaryConditionAdapter.h"
#include "StringUtil.hpp"
//#include "rapidjson/document.h"		// rapidjson's DOM-style API
//#include "rapidjson/filestream.h"

#include <fstream>
#include "yaml-cpp/yaml.h"

using namespace std;


void run(const char *istr)
{
   try
   {
      //// Prepare reader and input stream.
      //rapidjson::Reader reader;
      ////rapidjson::Document reader;
      //FILE* fp;
      //fp = fopen(istr, "r");
      //rapidjson::FileStream is(fp);

      //rapidjson::Document document;	// Default template parameter uses UTF8 and MemoryPoolAllocator.

      //if (document.ParseStream<0>(is).HasParseError())
      //{
      //   //UBLOG(logINFO,"JSON parcing is fail" );
      //   fprintf(stderr, "\nError(%u): %s\n", (unsigned)document.GetErrorOffset(), document.GetParseError());
      //   return;
      //}

      //fclose(fp);

      std::ifstream fin(istr);
      YAML::Parser parser(fin);
      YAML::Node doc;
      parser.GetNextDocument(doc);

      if(doc.FindValue("GeoObjects"))
      {
         const YAML::Node& geoObjects = doc["GeoObjects"];
         string id;
         for(unsigned i=0;i<geoObjects.size();i++)
         {
            geoObjects[i]["ID"] >> id;
            std::cout << id << "\n";
         }
      }

   }
   catch(YAML::ParserException& e) {
      std::cout << e.what() << "\n";
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
int main(int argc, char* argv[])
{
   if ( argv != NULL )
   {
      if (argc > 1)
      {
         run(argv[1]);
      }
      else
      {
         cout << "Input file must be set!: " <<  argv[0] << " <input file>" << endl << std::flush;
      }
   }

   return 0;
}


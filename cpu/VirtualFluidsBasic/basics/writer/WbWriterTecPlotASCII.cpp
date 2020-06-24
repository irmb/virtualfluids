#include <basics/writer/WbWriterTecPlotASCII.h>
#include <basics/utilities/UbLogger.h>

using namespace std;

/*===============================================================================*/
string WbWriterTecPlotASCII::writeOctsWithNodeData(const string& filename,vector<UbTupleFloat3 >& nodes, vector<UbTupleInt8 >& cells, vector<string >& datanames, vector<vector<double > >& nodedata)
{
   string tecplotfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterTecPlotASCII::writeOctsWithNodeData to "<<tecplotfilename<<" - start");

   ofstream out(tecplotfilename.c_str());
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(tecplotfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(tecplotfilename.c_str());}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+tecplotfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   out<<"TITLE = VirtualFluids OctGrid from "<<UbSystem::getTimeStamp()<<endl;

   out<<"VARIABLES = \"X\", \"Y\", \"Z\"";
   for(size_t d=0; d<datanames.size(); d++)   
      out<<", \""<<datanames[d]<<"\"";
   out<<endl;

   out<<"ZONE NODES="<<nofNodes<<", ELEMENTS="<<nofCells<<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<<endl;
   for(size_t n=0; n<nodes.size(); n++)   
   {
      UbTupleFloat3& coords = nodes[n];
      out<<val<1>(coords)<<" "
         <<val<2>(coords)<<" "
         <<val<3>(coords);
      for(size_t d=0; d<datanames.size(); d++)   
         out<<" "<<nodedata[d][n];
      out<<endl;
   }

   for(size_t c=0; c<cells.size(); c++)   
   {
      UbTupleInt8& cell = cells[c];
      out<<val<1>(cell)<<" "
         <<val<2>(cell)<<" "
         <<val<3>(cell)<<" "
         <<val<4>(cell)<<" "
         <<val<5>(cell)<<" "
         <<val<6>(cell)<<" "
         <<val<7>(cell)<<" "
         <<val<8>(cell)<<endl;
   }

   out.close();
   UBLOG(logDEBUG1,"WbWriterTecPlotASCII::writeOctsWithNodeData to "<<tecplotfilename<<" - end");

   return tecplotfilename;
}
/*===============================================================================*/
string WbWriterTecPlotASCII::writeOcts(const string& filename,vector< UbTupleFloat3 >& nodes, vector< UbTupleInt8 >& cells)
{
   vector<string > datanames;
   vector<vector<double > > nodedata;
   return writeOctsWithNodeData(filename,nodes,cells,datanames,nodedata);
}
/*===============================================================================*/

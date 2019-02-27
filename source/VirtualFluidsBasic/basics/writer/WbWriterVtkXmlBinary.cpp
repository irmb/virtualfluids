#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/utilities/UbLogger.h>
#include <cstring>

using namespace std;

/*===============================================================================*/
string WbWriterVtkXmlBinary::pvdEndTag ="   </Collection>\n</VTKFile>";
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeCollection(const string& filename, const vector<string>& filenames, const double& timeStep, const bool& sepGroups)
{
   string vtkfilename=filename+".pvd";
   ofstream out(vtkfilename.c_str());
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str());}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   string endian;
   if(UbSystem::isLittleEndian()) endian = "LittleEndian";
   else                           endian = "BigEndian";
   out<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\""<<endian<<"\" >"<<endl;
   out<<"   <Collection>"<<endl;
   
   int group = 0, part=0;
   for(size_t i=0; i<filenames.size(); i++)
   {
      out<<"       <DataSet timestep=\""<<timeStep<<"\" group=\""<<group<<"\" part=\""<<part<<"\" file=\""<<filenames[i]<<"\"/>"<<endl;
      if(sepGroups) group++;
      else          part++;
   }
   out<<pvdEndTag;
   out.close();

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::addFilesToCollection(const string& filename, const vector<string>& filenames, const double& timeStep, const bool& sepGroups)
{
   string vtkfilename=filename;
   fstream test(vtkfilename.c_str(), ios::in);
   if(!test)
   {
      test.clear();
      vtkfilename += ".pvd";
      test.open(vtkfilename.c_str(), ios::in);
      if(!test) return this->writeCollection(filename,filenames,timeStep,sepGroups);
   }
 
   fstream out(vtkfilename.c_str(), ios::in | ios::out);
   out.seekp(-(int)pvdEndTag.size()-1, ios_base::end);

   int group = 0;
   for(size_t i=0; i<filenames.size(); i++)
   {
      out<<"       <DataSet timestep=\""<<timeStep<<"\" group=\""<<group<<"\" part=\""<<i<<"\" file=\""<<filenames[i]<<"\"/>"<<endl;
      if(sepGroups) group++;
   }
   out<<pvdEndTag;

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeParallelFile(const string& filename,vector<string>& pieceSources, vector<string>& pointDataNames, vector<string>& cellDataNames)
{
   string vtkfilename=filename+".pvtu";
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeParallelFile to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str());
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str());}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   //VTK FILE
   out<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<"\n";
   out<<"  <PUnstructuredGrid GhostLevel=\"0\">"<<"\n";
   out<<"    <PPoints>\n"; 
   out<<"      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
   out<<"    </PPoints>\n";
   out<<"    <PPointData>\n";
   for(size_t s=0; s<pointDataNames.size(); s++)
      out<< "      <PDataArray type=\"Float32\" Name=\""<< pointDataNames[s] <<"\"/>\n";
   out<<"    </PPointData>\n";
   if (cellDataNames.size() > 0)
   {
      out<<"    <PCellData>\n";
      for(size_t s=0; s<cellDataNames.size(); s++)
         out<< "      <PDataArray type=\"Float32\" Name=\""<< cellDataNames[s] <<"\"/>\n";
      out<<"    </PCellData>\n";
   }

   for(size_t s=0; s<pieceSources.size(); s++)
      out<<"    <Piece Source=\""<<pieceSources[s]<<"\"/>\n";
   out<<"  </PUnstructuredGrid>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeParallelFile to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeLines(const string& filename,vector<UbTupleFloat3 >& nodes, vector<UbTupleInt2 >& lines)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeLines to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)lines.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 2 /*nodes per line */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per line */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of line */ * nofCells * sizeof(unsigned char);

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(lines[c]), sizeof(int) );
      out.write( (char*)&val<2>(lines[c]), sizeof(int) );
      
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 2 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 3;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeLines to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlBinary::writeLinesWithNodeData(const string& filename,vector<UbTupleFloat3 >& nodes, vector<UbTupleInt2 >& lines, std::vector< std::string >& datanames, std::vector< std::vector< double > >& nodedata)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeLinesWithNodeData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)lines.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 2 /*nodes per line  */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per line */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of line    */ * nofCells * sizeof(unsigned char);
   int bytesScalarData      = 1 /*scalar          */ * nofNodes * sizeof(float); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   //DATA SECTION
   out<<"         <PointData>\n";
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float32\" Name=\""<< datanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarData);
   }
   out<<"         </PointData>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(lines[c]), sizeof(int) );
      out.write( (char*)&val<2>(lines[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 3 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 5;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   //DATA SECTION
   //scalarData
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out.write((char*)&bytesScalarData,bytesPerByteVal);
      for(size_t d=0; d<nodedata[s].size(); ++d)
      {
         //loake kopie machen, da in nodedata "doubles" sind
         float tmp = (float)nodedata[s][d];
         out.write((char*)&tmp,sizeof(float));
      }
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeLinesWithNodeData to "<<vtkfilename<<" - end");

   return vtkfilename;

}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeTriangles(const string& filename,vector<UbTupleFloat3 >& nodes, vector<UbTupleInt3 >& triangles)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeTriangles to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }


   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)triangles.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3 - coord    */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 3 /*nodes per triangle  */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per triangle */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of triangle    */ * nofCells * sizeof(unsigned char);

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(triangles[c]), sizeof(int) );
      out.write( (char*)&val<2>(triangles[c]), sizeof(int) );
      out.write( (char*)&val<3>(triangles[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 3 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 5;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out<<flush;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeTriangles to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeTrianglesWithNodeData(const string& filename,vector< UbTupleFloat3 >& nodes, vector< UbTupleInt3 >& cells, vector< string >& datanames, vector< vector< double > >& nodedata)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeTrianglesWithNodeData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 3 /*nodes per tri   */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per tri  */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of tri     */ * nofCells * sizeof(unsigned char);
   int bytesScalarData      = 1 /*scalar          */ * nofNodes * sizeof(float); 
   
   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   //DATA SECTION
   out<<"         <PointData>\n";
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float32\" Name=\""<< datanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarData);
   }
   out<<"         </PointData>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 3 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 5;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   //DATA SECTION
   //scalarData
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out.write((char*)&bytesScalarData,bytesPerByteVal);
      for(size_t d=0; d<nodedata[s].size(); ++d)
      {
         //loake kopie machen, da in nodedata "doubles" sind
         float tmp = (float)nodedata[s][d];
         out.write((char*)&tmp,sizeof(float));
      }
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeTrianglesWithNodeData to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeQuads(const string& filename,vector< UbTupleFloat3 >& nodes, vector< UbTupleInt4 >& cells)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuads to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 4 /*nodes per quad  */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per quad */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of quad    */ * nofCells * sizeof(unsigned char);
  
   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<4>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 4 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 8;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out<<flush;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuads to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeQuadsWithNodeData(const string& filename,vector< UbTupleFloat3 >& nodes, vector< UbTupleInt4 >& cells, vector< string >& datanames, vector< vector< double > >& nodedata)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuadsWithNodeData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 4 /*nodes per quad  */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per quad */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of quad    */ * nofCells * sizeof(unsigned char);
   int bytesScalarData      = 1 /*scalar          */ * nofNodes * sizeof(float); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   //DATA SECTION
   out<<"         <PointData>\n";
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float64\" Name=\""<< datanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarData);
   }
   out<<"         </PointData>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<4>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 4 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 8;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   //DATA SECTION
   //scalarData
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out.write((char*)&bytesScalarData,bytesPerByteVal);
      for(size_t d=0; d<nodedata[s].size(); ++d)
      {
         //loake kopie machen, da in nodedata "doubles" sind
         float tmp = (float)nodedata[s][d];
         out.write((char*)&tmp,sizeof(float));
      }
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuadsWithNodeData to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeQuadsWithCellData(const string& filename,vector< UbTupleFloat3 >& nodes, vector< UbTupleInt4 >& cells, vector< string >& datanames, vector< vector< double > >& celldata)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuadsWithCellData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 4 /*nodes per quad  */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per quad */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of quad    */ * nofCells * sizeof(unsigned char);
   int bytesScalarData      = 1 /*scalar          */ * nofCells * sizeof(float); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   //DATA SECTION
   out<<"         <CellData>\n";
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float32\" Name=\""<< datanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarData);
   }
   out<<"         </CellData>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<4>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 4 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 8;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   //DATA SECTION
   //scalarData
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out.write((char*)&bytesScalarData,bytesPerByteVal);
      for(size_t d=0; d<celldata[s].size(); ++d)
      {
         //loake kopie machen, da in celldata "doubles" sind
         float tmp = (float)celldata[s][d];
         out.write((char*)&tmp,sizeof(float));
      }
   }

   out<<"\n</AppendedData>\n";

   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuadsWithCellData to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeQuadsWithNodeAndCellData(const string& filename,vector< UbTupleFloat3 >& nodes, vector< UbTupleInt4 >& cells, 
                                                                vector< string >& nodedatanames, vector< vector< double > >& nodedata, vector< string >& celldatanames,
                                                                vector< vector< double > >& celldata                                                                    )
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuadsWithNodeAndCellData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 4 /*nodes per quad  */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per quad */ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of quad    */ * nofCells * sizeof(unsigned char);
   int bytesScalarDataPoint = 1 /*scalar          */ * nofNodes * sizeof(float); 
   int bytesScalarDataCell  = 1 /*scalar          */ * nofCells * sizeof(float); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   // Point DATA SECTION
   out<<"         <PointData>\n";
   for(size_t s=0; s<nodedatanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float32\" Name=\""<< nodedatanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarDataPoint);
   }
   out<<"         </PointData>\n";


   // Cell DATA SECTION
   out<<"         <CellData>\n";
   for(size_t s=0; s<celldatanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float32\" Name=\""<< celldatanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarDataCell);
   }
   out<<"         </CellData>\n";
   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<4>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 4 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 8;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   //Point DATA SECTION
   //scalarData
   for(size_t s=0; s<nodedatanames.size(); ++s)
   {
      out.write((char*)&bytesScalarDataPoint,bytesPerByteVal);
      for(size_t d=0; d<nodedata[s].size(); ++d)
      {
         //loake kopie machen, da in nodedata "doubles" sind
         float tmp = (float)nodedata[s][d];
         out.write((char*)&tmp,sizeof(float));
      }
   }
   //Cell DATA SECTION
   //scalarData
   for(size_t s=0; s<celldatanames.size(); ++s)
   {
      out.write((char*)&bytesScalarDataCell,bytesPerByteVal);
      for(size_t d=0; d<celldata[s].size(); ++d)
      {
         //loake kopie machen, da in celldata "doubles" sind
         float tmp = (float)celldata[s][d];
         out.write((char*)&tmp,sizeof(float));
      }
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeQuadsWithNodeAndCellData to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeOctsWithCellData(const string& filename,vector<UbTupleFloat3 >& nodes, vector<UbTupleInt8 >& cells, vector<string >& datanames, vector<vector<double > >& celldata)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeOctsWithCellData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3      */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 8 /*nodes per oct */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per oct*/ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of oct   */ * nofCells * sizeof(unsigned char);
   int bytesScalarData      = 1 /*scalar        */ * nofCells * sizeof(float); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   //DATA SECTION
   out<<"         <CellData>\n";
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float32\" Name=\""<< datanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarData);
   }
   out<<"         </CellData>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<4>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
      out.write( (char*)&val<5>(cells[c]), sizeof(int) );
      out.write( (char*)&val<6>(cells[c]), sizeof(int) );
      out.write( (char*)&val<8>(cells[c]), sizeof(int) );
      out.write( (char*)&val<7>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 8 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 11;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   //DATA SECTION
   //scalarData
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out.write((char*)&bytesScalarData,bytesPerByteVal);
      for(size_t d=0; d<celldata[s].size(); ++d)
      {
         //loake kopie machen, da in celldata "doubles" sind
         float tmp = (float)celldata[s][d];
         out.write((char*)&tmp,sizeof(float));
      }
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeOctsWithCellData to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeOctsWithNodeData(const string& filename,vector<UbTupleFloat3 >& nodes, vector<UbTupleInt8 >& cells, vector<string >& datanames, vector<vector<double > >& nodedata)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeOctsWithNodeData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3      */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 8 /*nodes per oct */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per oct*/ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of oct   */ * nofCells * sizeof(unsigned char);
   int bytesScalarData      = 1 /*scalar        */ * nofNodes * sizeof(double); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"2.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   //DATA SECTION
   out<<"         <PointData>\n";
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float64\" Name=\""<< datanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarData);
   }
   out<<"         </PointData>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<4>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
      out.write( (char*)&val<5>(cells[c]), sizeof(int) );
      out.write( (char*)&val<6>(cells[c]), sizeof(int) );
      out.write( (char*)&val<8>(cells[c]), sizeof(int) );
      out.write( (char*)&val<7>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 8 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 11;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }

   //DATA SECTION
   //scalarData
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out.write((char*)&bytesScalarData,bytesPerByteVal);
      for(size_t d=0; d<nodedata[s].size(); ++d)
      {
         //loake kopie machen, da in nodedata "doubles" sind
         //float tmp = (float)nodedata[s][d];
         //out.write((char*)&tmp,sizeof(float));
         double tmp = nodedata[s][d];
         out.write((char*)&tmp,sizeof(double));
      }
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeOctsWithNodeData to "<<vtkfilename<<" - end");

   return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlBinary::writeOcts(const string& filename,vector< UbTupleFloat3 >& nodes, vector< UbTupleInt8 >& cells)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeOcts to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 
   int nofCells = (int)cells.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3      */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 8 /*nodes per oct */ * nofCells * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per oct*/ * nofCells * sizeof(int  );
   int bytesCellTypes       = 1 /*type of oct   */ * nofCells * sizeof(unsigned char);
   //int bytesScalarData      = 1 /*scalar        */ * nofNodes * sizeof(float); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofCells<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofCells; c++) 
   {
      out.write( (char*)&val<1>(cells[c]), sizeof(int) );
      out.write( (char*)&val<2>(cells[c]), sizeof(int) );
      out.write( (char*)&val<4>(cells[c]), sizeof(int) );
      out.write( (char*)&val<3>(cells[c]), sizeof(int) );
      out.write( (char*)&val<5>(cells[c]), sizeof(int) );
      out.write( (char*)&val<6>(cells[c]), sizeof(int) );
      out.write( (char*)&val<8>(cells[c]), sizeof(int) );
      out.write( (char*)&val<7>(cells[c]), sizeof(int) );
   }

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   int itmp;
   for(int c=1; c<=nofCells; c++)
   {
      itmp = 8 * c;    
      out.write( (char*)&itmp, sizeof(int) );
   }

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 11;
   for(int c=0; c<nofCells; c++)
   {
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeOcts to "<<vtkfilename<<" - end");

   return vtkfilename;
}
std::string WbWriterVtkXmlBinary::writeNodes(const std::string& filename,std::vector< UbTupleFloat3 >& nodes)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeNodes to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 1 /*nodes per cell  */ * nofNodes * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per cell */ * nofNodes * sizeof(int  );
   int bytesCellTypes       = 1 /*type of line    */ * nofNodes * sizeof(unsigned char);

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofNodes<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofNodes; c++) 
      out.write( (char*)&c, sizeof(int) );

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   for(int c=1; c<=nofNodes; c++)
      out.write( (char*)&c, sizeof(int) );

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 1;
   for(int c=0; c<nofNodes; c++)
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );

   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeNodes to "<<vtkfilename<<" - end");

   return vtkfilename;
}
std::string WbWriterVtkXmlBinary::writeNodesWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata)
{
   string vtkfilename = filename+getFileExtension();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeNodesWithNodeData to "<<vtkfilename<<" - start");

   ofstream out(vtkfilename.c_str(),ios::out | ios::binary);
   if(!out)
   { 
      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
      string path = UbSystem::getPathFromString(vtkfilename);
      if(path.size()>0){ UbSystem::makeDirectory(path); out.open(vtkfilename.c_str(),ios::out | ios::binary);}
      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
   }

   int nofNodes = (int)nodes.size(); 

   int bytesPerByteVal      = 4; //==sizeof(int)
   int bytesPoints          = 3 /*x1/x2/x3       */ * nofNodes * sizeof(float);
   int bytesCellConnectivty = 1 /*nodes per cell */ * nofNodes * sizeof(int  );
   int bytesCellOffsets     = 1 /*offset per cell*/ * nofNodes * sizeof(int  );
   int bytesCellTypes       = 1 /*type of oct    */ * nofNodes * sizeof(unsigned char);
   int bytesScalarData      = 1 /*scalar         */ * nofNodes * sizeof(double); 

   int offset = 0;
   //VTK FILE
   out<<"<?xml version=\"1.0\"?>\n";
   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
   out<<"   <UnstructuredGrid>"<<"\n";
   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\" NumberOfCells=\""<<nofNodes<<"\">\n";

   //POINTS SECTION
   out<<"         <Points>\n"; 
   out<<"            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<< offset <<"\"  />\n";
   out<<"         </Points>\n";
   offset += (bytesPerByteVal + bytesPoints);

   //CELLS SECTION
   out<<"         <Cells>\n";
   out<<"            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellConnectivty); 
   out<<"            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<< offset <<"\" />\n";
   offset += (bytesPerByteVal + bytesCellOffsets);
   out<<"            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<< offset <<"\" />\n ";
   offset += (bytesPerByteVal + bytesCellTypes);
   out<<"         </Cells>\n";

   //DATA SECTION
   out<<"         <PointData>\n";
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out<< "            <DataArray type=\"Float64\" Name=\""<< datanames[s] <<"\" format=\"appended\" offset=\""<< offset <<"\" /> \n";
      offset += (bytesPerByteVal + bytesScalarData);
   }
   out<<"         </PointData>\n";

   out<<"      </Piece>\n";
   out<<"   </UnstructuredGrid>\n";

   // AppendedData SECTION
   out<<"   <AppendedData encoding=\"raw\">\n";
   out<<"_";

   //POINTS SECTION
   out.write((char*)&bytesPoints,bytesPerByteVal);
   for(int n=0; n<nofNodes; n++)
   {
      out.write((char*)&val<1>(nodes[n]),sizeof(float));
      out.write((char*)&val<2>(nodes[n]),sizeof(float));
      out.write((char*)&val<3>(nodes[n]),sizeof(float));
   }

   //CELLS SECTION
   //cellConnectivity
   out.write( (char*)&bytesCellConnectivty, bytesPerByteVal );  
   for(int c=0; c<nofNodes; c++) 
      out.write( (char*)&c, sizeof(int) );

   //cellOffsets
   out.write( (char*)&bytesCellOffsets, bytesPerByteVal );
   for(int c=1; c<=nofNodes; c++)
      out.write( (char*)&c, sizeof(int) );

   //cellTypes
   out.write( (char*)&bytesCellTypes, bytesPerByteVal );
   unsigned char vtkCellType = 1;
   for(int c=0; c<nofNodes; c++)
      out.write( (char*)&vtkCellType, sizeof(unsigned char) );
   
   //DATA SECTION
   //scalarData
   for(size_t s=0; s<datanames.size(); ++s)
   {
      out.write((char*)&bytesScalarData,bytesPerByteVal);
      for(size_t d=0; d<nodedata[s].size(); ++d)
      {
         //loake kopie machen, da in nodedata "doubles" sind
         //float tmp = (float)nodedata[s][d];
         //out.write((char*)&tmp,sizeof(float));
         double tmp = nodedata[s][d];
         out.write((char*)&tmp, sizeof(double));
      }
   }
   out<<"\n</AppendedData>\n";
   out<<"</VTKFile>";
   out<<endl;
   out.close();
   UBLOG(logDEBUG1,"WbWriterVtkXmlBinary::writeNodesWithNodeData to "<<vtkfilename<<" - end");

   return vtkfilename;

}

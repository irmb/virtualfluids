#ifdef CAB_ZLIB
   #include <basics/writer/WbWriterBOBJ.h>
   #include <basics/utilities/UbLogger.h>
   #include <cstring>

   #include <zlib.h>


   using namespace std;
   /*===============================================================================*/
   std::string WbWriterBOBJ::writeTriangles(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt3 >& triangles)
   {
      string bobjFilename=filename+getFileExtension();
      UBLOG(logDEBUG1,"WbWriterBOBJ::writeTriangles to "<<bobjFilename<<" - start");

      gzFile gzf = gzopen( bobjFilename.c_str(), "wb1" );
      
      size_t nofNodes     = nodes.size(); 
      size_t nofTriangles = triangles.size(); 

      //write to file
      size_t numVerts;
      //double v[3];
      if(sizeof(numVerts)!=4) { throw UbException(UB_EXARGS,"danger..."); }
      numVerts = nofNodes;
      gzwrite(gzf, &numVerts, sizeof(numVerts));

      for(size_t k=0; k<nofNodes; k++) {
         float vertp = val<1>(nodes[k]);
         gzwrite(gzf, &vertp, sizeof(vertp));
         vertp       = val<2>(nodes[k]);
         gzwrite(gzf, &vertp, sizeof(vertp));
         vertp       = val<3>(nodes[k]);
         gzwrite(gzf, &vertp, sizeof(vertp));
      }

      //NORMAL VECTOR
      //double n[3];
      gzwrite(gzf, &numVerts, sizeof(numVerts));
      for(size_t k=0; k<nofNodes; k++) {
         //poly->GetPointData()->GetNormals()->GetTuple(k, n);
         float normp = 0.0;//n[0];
         gzwrite(gzf, &normp, sizeof(normp));
         normp = 0.0;//n[1];
         gzwrite(gzf, &normp, sizeof(normp));
         normp = 0.0;//n[2];
         gzwrite(gzf, &normp, sizeof(normp));
      }

      //vtkIdType npts = 3;
      //vtkIdType* pts;
      size_t numTris = nofTriangles;
      gzwrite(gzf, &numTris, sizeof(numTris));
      for(size_t k=0; k<nofTriangles/*(size_t)poly->GetNumberOfPolys()*/; k++) {
         //poly->GetPolys()->GetNextCell(npts, pts);
         //int triIndex = *pts;
         //gzwrite(gzf, &triIndex, sizeof(triIndex)); 
         //triIndex = *(pts+1);
         //gzwrite(gzf, &triIndex, sizeof(triIndex)); 
         //triIndex = *(pts+2);
         //gzwrite(gzf, &triIndex, sizeof(triIndex));
         //poly->GetPolys()->GetNextCell(npts, pts);
         int triIndex = val<1>(triangles[k]);//*pts;
         gzwrite(gzf, &triIndex, sizeof(triIndex)); 
         triIndex     = val<2>(triangles[k]);//*(pts+1);
         gzwrite(gzf, &triIndex, sizeof(triIndex)); 
         triIndex     = val<3>(triangles[k]);//*(pts+2);
         gzwrite(gzf, &triIndex, sizeof(triIndex));
      }

      gzclose( gzf );

      UBLOG(logDEBUG1,"WbWriterBOBJ::writeTriangles to "<<bobjFilename<<" - end");

      return bobjFilename;
   }
   /*===============================================================================*/

#endif //CAB_ZLIB

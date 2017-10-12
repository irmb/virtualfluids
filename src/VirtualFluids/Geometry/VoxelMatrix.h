//template <class T>
//struct OriginVoxelMatrix
//{
//   CbArray3D<T> matrix;
//};
////////////////////////////////////////////////////////////////////////////
//template <class T>
//void GbVoxelMatrix3D::readOriginMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian)
//{
//   using namespace std;
//   ifstream in(filename.c_str(), ios::binary);
//   if (!in) throw UbException(UB_EXARGS, "could not open file "+filename);
//
//   in.seekg(0, ios::end);     //Ende springen
//   fstream::off_type length = in.tellg(); //Position abfragen
//   in.seekg(0, ios::beg);    //An den Anfang springen 
//
//   unsigned long long nofn = (unsigned long long)nodesX1*(unsigned long long)nodesX2*(unsigned long long)nodesX3*(unsigned long long)sizeof(T);
//   if (nofn!=(unsigned long long)length)
//   {
//      throw UbException(UB_EXARGS, "number of nodes("+UbSystem::toString(nofn)+") doesn't match file size("+UbSystem::toString((long)length)+")");
//   }
//
//   GbVoxelMatrix3D::oVoxelMatrix3D<T> = CbArray3D<T>(nodesX1, nodesX2, nodesX3);
//
//   T val;
//   for (int x3 = 0; x3<nodesX3; x3++)
//      for (int x2 = 0; x2<nodesX2; x2++)
//         for (int x1 = 0; x1<nodesX1; x1++)
//         {
//            in.read((char*)&val, sizeof(T));
//            if (endian==BigEndian)
//               UbSystem::swapByteOrder((unsigned char*)(&(val)), sizeof(T));
//            oVoxelMatrix3D<T>(x1, x2, x3) = val;
//         }
//}
////////////////////////////////////////////////////////////////////////////
//template <class T>
//void GbVoxelMatrix3D::writeOriginMatrixFromRawFile(std::string filename)
//{
//   using namespace std;
//
//   string fn = filename;
//
//   FILE *file;
//   file = fopen(fn.c_str(), "wb");
//
//   if (file==NULL)
//   {
//      std::string pathf = UbSystem::getPathFromString(fn);
//      if (fn.size()>0) { UbSystem::makeDirectory(pathf); file = fopen(fn.c_str(), "wb"); }
//      if (file==NULL) throw UbException(UB_EXARGS, "can not open "+fn);
//   }
//
//   fwrite(GbVoxelMatrix3D::oVoxelMatrix3D<T>.getStartAdressOfSortedArray(0, 0, 0), sizeof(T), GbVoxelMatrix3D::oVoxelMatrix3D<T>.getDataVector().size(), file);
//   fclose(file);
//}

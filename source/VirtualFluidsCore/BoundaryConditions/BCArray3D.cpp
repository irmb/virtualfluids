#include "BCArray3D.h"


const int BCArray3D::SOLID = -1;
const int BCArray3D::FLUID = -2;
const int BCArray3D::INTERFACECF = -3;
const int BCArray3D::INTERFACEFC = -4;
const int BCArray3D::UNDEFINED = -5;

//////////////////////////////////////////////////////////////////////////
BCArray3D::BCArray3D() {}
//////////////////////////////////////////////////////////////////////////
BCArray3D::BCArray3D(std::size_t nx1, std::size_t nx2, std::size_t nx3)
{
   bcindexmatrix.resize(nx1, nx2, nx3, UNDEFINED);
}
//////////////////////////////////////////////////////////////////////////
BCArray3D::BCArray3D(std::size_t nx1, std::size_t nx2, std::size_t nx3, int val)
{
   bcindexmatrix.resize(nx1, nx2, nx3, val);
}
//////////////////////////////////////////////////////////////////////////
BCArray3D::~BCArray3D() {}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::resize(std::size_t nx1, std::size_t nx2, std::size_t nx3)
{
   bcindexmatrix.resize(nx1, nx2, nx3);
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::resize(std::size_t nx1, std::size_t nx2, std::size_t nx3, int val)
{
   bcindexmatrix.resize(nx1, nx2, nx3, val);
}
//////////////////////////////////////////////////////////////////////////
bool BCArray3D::validIndices(std::size_t x1, std::size_t x2, std::size_t x3)  const
{
   if (x1 < 0 || x1 >= this->bcindexmatrix.getNX1()) return false;
   if (x2 < 0 || x2 >= this->bcindexmatrix.getNX2()) return false;
   if (x3 < 0 || x3 >= this->bcindexmatrix.getNX3()) return false;
   return true;
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setBC(std::size_t x1, std::size_t x2, std::size_t x3, BCClassPtr const& bc)
{
   if (this->hasBC(x1, x2, x3))
   {
      if (this->getBC(x1, x2, x3) == bc) return;
      else                            this->deleteBC(x1, x2, x3);
   }

   //wenn keine frei gewordene BCs vorhanden
   if (indexContainer.empty())
   {
      bcvector.push_back(bc);
      bcindexmatrix(x1, x2, x3) = (int)bcvector.size() - 1;
   }
   else
   {
      int index = indexContainer.back();
      bcindexmatrix(x1, x2, x3) = index;
      bcvector[index] = bc;
      indexContainer.pop_back();
   }
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setSolid(std::size_t x1, std::size_t x2, std::size_t x3)
{
   if (this->hasBC(x1, x2, x3)) this->deleteBC(x1, x2, x3);
   bcindexmatrix(x1, x2, x3) = SOLID;
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setFluid(std::size_t x1, std::size_t x2, std::size_t x3)
{
   if (this->hasBC(x1, x2, x3)) this->deleteBC(x1, x2, x3);
   bcindexmatrix(x1, x2, x3) = FLUID;
}
//////////////////////////////////////////////////////////////////////////
void BCArray3D::setUndefined(std::size_t x1, std::size_t x2, std::size_t x3)
{
   if (this->hasBC(x1, x2, x3)) this->deleteBC(x1, x2, x3);
   bcindexmatrix(x1, x2, x3) = UNDEFINED;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfSolidEntries() const
{
   const std::vector<int>& data = bcindexmatrix.getDataVector();
   std::size_t counter = 0;
   for (std::size_t i = 0; i < data.size(); i++)
      if (data[i] == SOLID) counter++;
   return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfFluidEntries() const
{
   const std::vector<int>& data = bcindexmatrix.getDataVector();
   std::size_t counter = 0;
   for (std::size_t i = 0; i < data.size(); i++)
   {
      int tmp = data[i];
      if (tmp == FLUID || tmp >= 0) counter++;
   }
   return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfFluidWithoutBCEntries() const
{
   const std::vector<int>& data = bcindexmatrix.getDataVector();
   std::size_t counter = 0;
   for (std::size_t i = 0; i < data.size(); i++)
      if (data[i] == FLUID) counter++;
   return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfBCEntries() const
{
   const std::vector<int>& data = bcindexmatrix.getDataVector();
   std::size_t counter = 0;
   for (std::size_t i = 0; i < data.size(); i++)
      if (data[i] >= 0) counter++;
   return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getNumberOfUndefinedEntries() const
{
   const std::vector<int>& data = bcindexmatrix.getDataVector();
   std::size_t counter = 0;
   for (std::size_t i = 0; i < data.size(); i++)
      if (data[i] == UNDEFINED) counter++;
   return counter;
}
//////////////////////////////////////////////////////////////////////////
std::size_t BCArray3D::getBCVectorSize() const
{
   return this->bcvector.size();
}
//////////////////////////////////////////////////////////////////////////
std::string BCArray3D::toString() const
{
   std::size_t solidCounter = 0;
   std::size_t fluidCounter = 0;
   std::size_t bcCounter = 0;
   std::size_t undefCounter = 0;

   for (int x1 = 0; x1 < bcindexmatrix.getNX1(); x1++)
   {
      for (int x2 = 0; x2 < bcindexmatrix.getNX2(); x2++)
      {
         for (int x3 = 0; x3 < bcindexmatrix.getNX3(); x3++)
         {
            if (bcindexmatrix(x1, x2, x3) >= 0) bcCounter++;
            else if (bcindexmatrix(x1, x2, x3) == FLUID) fluidCounter++;
            else if (bcindexmatrix(x1, x2, x3) == SOLID) solidCounter++;
            else if (bcindexmatrix(x1, x2, x3) == UNDEFINED) undefCounter++;
            else throw UbException(UB_EXARGS, "invalid matrixEntry");
         }
      }
   }

   std::size_t unrefEntriesInBcVector = 0;
   for (std::size_t i = 0; i < bcvector.size(); i++) if (!bcvector[i]) unrefEntriesInBcVector++;

   std::stringstream text;
   text << "BCArray<" << typeid(BCClassPtr).name() << "," << typeid(int).name() << ">";
   text << "[ entries: " << bcindexmatrix.getNX1() << "x" << bcindexmatrix.getNX2();
   text << "x" << bcindexmatrix.getNX3() << "=";
   text << bcindexmatrix.getNX1()*bcindexmatrix.getNX2()*bcindexmatrix.getNX3() << " ]:\n";
   text << " - #fluid entries : " << fluidCounter << std::endl;
   text << " - #bc    entries : " << bcCounter << std::endl;
   text << " - #solid entries : " << solidCounter << std::endl;
   text << " - #undef entries : " << undefCounter << std::endl;
   text << " - bcvector-entries      : " << bcvector.size() << " (empty ones: " << unrefEntriesInBcVector << ")\n";
   text << " - indexContainer-entries: " << indexContainer.size() << std::endl;

   return text.str();
}
//////////////////////////////////////////////////////////////////////////
std::vector< int >& BCArray3D::getBcindexmatrixDataVector()
{
   return bcindexmatrix.getDataVector();
}


//////////////////////////////////////////////////////////////////////////
void BCArray3D::deleteBCAndSetType(std::size_t x1, std::size_t x2, std::size_t x3, int type)
   {
      this->deleteBC(x1, x2, x3);

      //matrix neuen Typ zuweisen
      bcindexmatrix(x1, x2, x3) = type;
   }
//////////////////////////////////////////////////////////////////////////
void BCArray3D::deleteBC(std::size_t x1, std::size_t x2, std::size_t x3)
   {
      //ueberpruefen, ob ueberhaupt BC vorhanden
      int index = bcindexmatrix(x1, x2, x3);
      if (index < 0) return;

      //frei gewordenen Index in den Indexcontainer schieben
      indexContainer.push_back(index);

      //element "loeschen"
      bcvector[index] = BCClassPtr();
   }
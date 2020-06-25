#include "kFullReader.h"

#include "VirtualFluidsBasics/basics/utilities/UbFileInputASCII.h"

void kFullReader::readFileForAlloc(const std::string fileName, Parameter* para)
{
   UbFileInputASCII in(fileName);

   int maxlevel = in.readInteger();
   in.readLine();
   
   for (int level = 0; level <= maxlevel; level++)
   {
	   int nodeNumberX1 = in.readInteger() - (2 * level);
	   int nodeNumberX2 = in.readInteger() - (2 * level);
	   int nodeNumberX3 = in.readInteger() - (2 * level);
	   para->getParH(level)->gridNX = nodeNumberX1;
	   para->getParH(level)->gridNY = nodeNumberX2;
	   para->getParH(level)->gridNZ = nodeNumberX3;
	   para->getParD(level)->gridNX = para->getParH(level)->gridNX;
	   para->getParD(level)->gridNY = para->getParH(level)->gridNY;
	   para->getParD(level)->gridNZ = para->getParH(level)->gridNZ;

	   //weitere Werte setzen
	   para->getParH(level)->nx                    = para->getParH(level)->gridNX + 2 * STARTOFFX;
	   para->getParH(level)->ny                    = para->getParH(level)->gridNY + 2 * STARTOFFY;
	   para->getParH(level)->nz                    = para->getParH(level)->gridNZ + 2 * STARTOFFZ;
	   para->getParH(level)->size_Mat              = para->getParH(level)->nx * para->getParH(level)->ny * para->getParH(level)->nz;
	   para->getParH(level)->mem_size_int          = sizeof(unsigned int) * para->getParH(level)->size_Mat;
	   //
	   para->getParD(level)->nx                    = para->getParD(level)->gridNX + 2 * STARTOFFX;
	   para->getParD(level)->ny                    = para->getParD(level)->gridNY + 2 * STARTOFFY;
	   para->getParD(level)->nz                    = para->getParD(level)->gridNZ + 2 * STARTOFFZ;
	   para->getParD(level)->size_Mat              = para->getParD(level)->nx * para->getParD(level)->ny * para->getParD(level)->nz;
	   para->getParD(level)->mem_size_int          = sizeof(unsigned int) * para->getParD(level)->size_Mat;


	   in.readLine();
	   for(unsigned int ix3=STARTOFFZ - level; ix3<para->getParH(level)->gridNZ + STARTOFFZ + level ; ix3++)
	   {
		   for(unsigned int ix2=STARTOFFY - level; ix2<para->getParH(level)->gridNY + STARTOFFY + level ; ix2++)
		   {
			   for(unsigned int ix1=STARTOFFX - level; ix1<para->getParH(level)->gridNX + STARTOFFX + level ; ix1++)
			   {
				   in.readInteger();
			   }
		   }
	   }
	   in.readLine();
   } //end levelloop
}

void kFullReader::readFile(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);

	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level <= maxlevel; level++)
	{
		in.readInteger();
		in.readInteger();
		in.readInteger();

		in.readLine();
		for(unsigned int ix3=STARTOFFZ - level; ix3<para->getParH(level)->gridNZ + STARTOFFZ + level; ix3++)
		{
			for(unsigned int ix2=STARTOFFY - level; ix2<para->getParH(level)->gridNY + STARTOFFY + level; ix2++)
			{
				for(unsigned int ix1=STARTOFFX - level; ix1<para->getParH(level)->gridNX + STARTOFFX + level; ix1++)
				{
					unsigned int m = para->getParH(level)->nx*(para->getParH(level)->ny*ix3 + ix2) + ix1;
					para->getParH(level)->k[m] =  in.readInteger();//+1;//!!achtung + 1
				}
			}
		}
		in.readLine();
	} //end levelloop
}


//GEO-FULL
void kFullReader::readGeoFull(const std::string fileName, Parameter* para)
{
	UbFileInputASCII in(fileName);
	int test = 0;
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level <= maxlevel; level++)
	{
		in.readInteger();
		in.readInteger();
		in.readInteger();

		in.readLine();
		for(unsigned int ix3=STARTOFFZ - level; ix3<para->getParH(level)->gridNZ + STARTOFFZ + level; ix3++)
		{
			for(unsigned int ix2=STARTOFFY - level; ix2<para->getParH(level)->gridNY + STARTOFFY + level; ix2++)
			{
				for(unsigned int ix1=STARTOFFX - level; ix1<para->getParH(level)->gridNX + STARTOFFX + level; ix1++)
				{
					unsigned int m = para->getParH(level)->nx*(para->getParH(level)->ny*ix3 + ix2) + ix1;
					test = in.readInteger();
					if (test == 7)
					{
						test = 1;
					}
					else if (test != 1)//???
					{
						test = 16;
					}
					para->getParH(level)->geo[m] = test;
					//para->getParH(level)->geo[m] =  in.readInteger();//+1;//!!achtung + 1
				}
			}
		}
		in.readLine();
	} //end levelloop
}
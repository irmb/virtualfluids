#include "PositionReader.h"

#include "Parameter/Parameter.h"

#include <basics/utilities/UbFileInputASCII.h>

using namespace vf::lbm::dir;

//////////////////////////////////////////////////////////////////////////
void PositionReader::readMeasurePoints( Parameter* para ) 
{
    UbFileInputASCII in(para->getmeasurePoints());
    int numberOfAllNodes = in.readInteger();
    in.readLine();
    int tempLevel;
    MeasurePoints tempMP;
    //printf("done, init the values...\n");
    for (int u = 0; u < numberOfAllNodes; u++)
    {
        tempMP.name = in.readString();         
        //printf("done, read the name...\n");
        tempMP.k = in.readInteger();
        //printf("done, read k...\n");
        tempLevel = in.readInteger();
        //printf("done, read level...\n");
        in.readLine();
        //printf("done, read the values...\n");
        para->getParH(tempLevel)->MP.push_back(tempMP);
        //printf("done, put it into a vector...\n");
    }
}

#ifndef FILE_WRITING_INFORMATION_STRUCT_H
#define FILE_WRITING_INFORMATION_STRUCT_H

#include <string>

struct VectorWriterInformationStruct
{
    unsigned int startTimeVectorWriter;
    bool writeVTKFiles;
    unsigned int startTimeVTKDataWriter;
};

#endif 
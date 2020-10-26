#ifndef VIRTUALFLUIDSPYTHONBINDINGS_WRITERCONFIG_H
#define VIRTUALFLUIDSPYTHONBINDINGS_WRITERCONFIG_H

#include <string>
#include <basics/writer/WbWriter.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

enum WriterType {
    ASCII, BINARY
};

struct WriterConfig {
    WriterType writerType{};
    std::string outputPath{"./output"};

    WbWriter *getWriter()
    {
        if (writerType == ASCII) return WbWriterVtkXmlASCII::getInstance();
        if (writerType == BINARY) return WbWriterVtkXmlBinary::getInstance();
        return nullptr;
    }
};

#endif //VIRTUALFLUIDSPYTHONBINDINGS_WRITERCONFIG_H

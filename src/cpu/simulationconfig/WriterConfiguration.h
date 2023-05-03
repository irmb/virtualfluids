#ifndef VIRTUALFLUIDSPYTHONBINDINGS_WRITERCONFIG_H
#define VIRTUALFLUIDSPYTHONBINDINGS_WRITERCONFIG_H

#include <string>
#include <basics/writer/WbWriter.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

enum OutputFormat {
    ASCII, BINARY
};

struct WriterConfiguration {
    OutputFormat outputFormat{};
    std::string outputPath{"./output"};

    WbWriter *getWriter() const 
    {
        if (outputFormat == ASCII) return WbWriterVtkXmlASCII::getInstance();
        if (outputFormat == BINARY) return WbWriterVtkXmlBinary::getInstance();
        return nullptr;
    }
};

#endif //VIRTUALFLUIDSPYTHONBINDINGS_WRITERCONFIG_H

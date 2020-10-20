#ifndef WBWRITERX3D_H
#define WBWRITERX3D_H

#include <string>

#include <basics/writer/WbWriter.h>

class WbWriterX3D : public WbWriter
{
public:
    static WbWriterX3D *getInstance()
    {
        static WbWriterX3D instance;
        return &instance;
    }

    WbWriterX3D(const WbWriterX3D &) = delete;
    const WbWriterX3D &operator=(const WbWriterX3D &) = delete;

private:
    WbWriterX3D() : WbWriter()
    {
        if (sizeof(unsigned char) != 1)
            throw UbException(UB_EXARGS, "error char  type mismatch");
        if (sizeof(int) != 4)
            throw UbException(UB_EXARGS, "error int   type mismatch");
        if (sizeof(float) != 4)
            throw UbException(UB_EXARGS, "error float type mismatch");
    }

    static std::string pvdEndTag;

public:
    std::string getFileExtension() override { return "ascii.X3D"; }

    std::string writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                               std::vector<UbTupleInt3> &triangles) override;
};

#endif // WBWRITERX3D_H

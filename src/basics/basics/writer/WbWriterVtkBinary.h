//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef WBWRITERVTKBINARY_H
#define WBWRITERVTKBINARY_H

#include <basics/writer/WbWriter.h>

class WbWriterVtkBinary : public WbWriter
{
public:
    static WbWriterVtkBinary *getInstance()
    {
        static WbWriterVtkBinary instance;
        return &instance;
    }

    WbWriterVtkBinary(const WbWriterVtkBinary &) = delete;
    const WbWriterVtkBinary &operator=(const WbWriterVtkBinary &) = delete;

private:
    WbWriterVtkBinary() = default;

public:
    std::string getFileExtension() override { return ".bin.vtk"; }

    //////////////////////////////////////////////////////////////////////////
    // lines
    //     0 ---- 1
    // nodenumbering must start with 0!
    std::string writeLines(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt2> &lines) override;

    //////////////////////////////////////////////////////////////////////////
    // triangles
    // cell numbering:
    //                    2
    //
    //                  0---1
    // nodenumbering must start with 0!
    std::string writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                               std::vector<UbTupleInt3> &cells) override;
    std::string writeTrianglesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                           std::vector<UbTupleInt3> &cells, std::vector<std::string> &datanames,
                                           std::vector<std::vector<double>> &nodedata) override;

    //////////////////////////////////////////////////////////////////////////
    // 2D
    // cell numbering:
    //                  3---2
    //                  |   |
    //                  0---1
    // nodenumbering must start with 0!
    std::string writeQuads(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt4> &cells) override;
    std::string writeQuadsWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                       std::vector<UbTupleInt4> &cells, std::vector<std::string> &datanames,
                                       std::vector<std::vector<double>> &nodedata) override;
    std::string writeQuadsWithCellData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                       std::vector<UbTupleInt4> &cells, std::vector<std::string> &datanames,
                                       std::vector<std::vector<double>> &celldata) override;
    std::string writeQuadsWithNodeAndCellData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                              std::vector<UbTupleInt4> &cells, std::vector<std::string> &nodedatanames,
                                              std::vector<std::vector<double>> &nodedata,
                                              std::vector<std::string> &celldatanames,
                                              std::vector<std::vector<double>> &celldata) override;

    //////////////////////////////////////////////////////////////////////////
    // octs
    //     7 ---- 6
    //    /|     /|
    //   4 +--- 5 |
    //   | |    | |
    //   | 3 ---+ 2
    //   |/     |/
    //   0 ---- 1
    std::string writeOcts(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                          std::vector<UbTupleInt8> &cells) override;
    std::string writeOctsWithCellData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                      std::vector<UbTupleInt8> &cells, std::vector<std::string> &datanames,
                                      std::vector<std::vector<double>> &celldata) override;
    std::string writeOctsWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                      std::vector<UbTupleUInt8> &cells, std::vector<std::string> &datanames,
                                      std::vector<std::vector<double>> &nodedata) override;
};

#endif // WBWRITERVTKBINARY_H

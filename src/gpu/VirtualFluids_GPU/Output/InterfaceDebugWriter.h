#ifndef INTERFACEDEBUG_HPP
#define INTERFACEDEBUG_HPP

class Parameter;
namespace InterfaceDebugWriter
{

void writeInterfaceLinesDebugCF(Parameter *para);
void writeInterfaceLinesDebugFC(Parameter *para);

void writeInterfaceLinesDebugCFCneighbor(Parameter *para);
void writeInterfaceLinesDebugCFFneighbor(Parameter *para);
void writeInterfaceLinesDebugFCCneighbor(Parameter *para);
void writeInterfaceLinesDebugFCFneighbor(Parameter *para);

void writeInterfaceLinesDebugOff(Parameter *para);
void writeInterfacePointsDebugCFC(Parameter *para);

void writeBcPointsDebug(Parameter *para);
void writePressPointsDebug(Parameter *para);
void writePressNeighborPointsDebug(Parameter *para);

void writeNeighborXPointsDebug(Parameter *para);
void writeNeighborXLinesDebug(Parameter *para);
void writeNeighborYPointsDebug(Parameter *para);
void writeNeighborYLinesDebug(Parameter *para);
void writeNeighborZPointsDebug(Parameter *para);
void writeNeighborZLinesDebug(Parameter *para);

void writeInterfaceCellsDebugCFC(Parameter *para);
void writeInterfaceCellsDebugCFF(Parameter *para);

void writeInterfaceFCC_Send(Parameter *para, int processID = 0);
void writeInterfaceCFC_Recv(Parameter *para, int processID = 0);

void writeSendNodesStream(Parameter *para, int processID = 0);
void writeRecvNodesStream(Parameter *para, int processID = 0);
} // namespace InterfaceDebugWriter

#endif
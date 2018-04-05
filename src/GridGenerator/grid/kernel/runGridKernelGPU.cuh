#ifndef runGridKernelGPU_H
#define runGridKernelGPU_H


class Grid;
class GridImp;
struct TriangularMesh;
class LaunchParameter;

float runKernelInitalUniformGrid3d(const LaunchParameter& para, GridImp &grid);
float runKernelToMesh(const LaunchParameter& para, GridImp &grid, const TriangularMesh &geom);
float runKernelToMarkNodesToDeleteOutsideOfGeometry(const LaunchParameter& para, GridImp &grid);

float runKernelToFindGridInterface(const LaunchParameter& para, GridImp &grid, GridImp &finerGrid);
float runKernelToFindNeighborsNewIndices(const LaunchParameter& para, GridImp &grid);
float runKernelToFindGridInterfaceNewIndices(const LaunchParameter& para, GridImp &grid);


float runKernelSetOverlapNodesToInvalid(const LaunchParameter& para, GridImp &grid, GridImp &finerGrid);
float runKernelSetToInvalid(const LaunchParameter& para, GridImp &grid);
float runKernelFindIndices(const LaunchParameter& para, GridImp &grid);

#endif

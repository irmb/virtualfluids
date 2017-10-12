#ifndef runGridKernelGPU_H
#define runGridKernelGPU_H


struct Grid;
struct Geometry;
class LaunchParameter;

float runKernelInitalUniformGrid3d(const LaunchParameter& para, Grid &grid);
float runKernelToMesh(const LaunchParameter& para, Grid &grid, const Geometry &geom);
float runKernelToMarkNodesToDeleteOutsideOfGeometry(const LaunchParameter& para, Grid &grid);

float runKernelSetToInvalid(const LaunchParameter& para, Grid &grid);
float runKernelFindIndices(const LaunchParameter& para, Grid &grid);

#endif

#ifndef geometry_helper_H
#define geometry_helper_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>


__host__ __device__ __inline__ void rotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY, real &originX, real &originY)
{
    real distX = posX - originX;
    real distY = posY - originY;

    newPosX = distX*cos(angle) - distY*sin(angle);
    newPosY = distX*sin(angle) + distY*cos(angle);  

    newPosX += originX;
    newPosY += originY;
}

__host__ __device__ __inline__ void invRotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY, real &originX, real &originY)
{
    real distX = posX - originX;
    real distY = posY - originY;

    newPosX =  distX*cos(angle) + distY*sin(angle);
    newPosY = -distX*sin(angle) + distX*cos(angle);  

    newPosX += originX;
    newPosY += originY;
}

__host__ __device__ __inline__ void rotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real distX = posX - originX;
    real distY = posY - originY;
    real distZ = posZ - originZ;

    newPosX = posX;
    newPosY = distY*cos(angle) - distZ*sin(angle);
    newPosZ = distY*sin(angle) + distZ*cos(angle);

    newPosX += originX;
    newPosY += originY;
    newPosZ += originZ;
}

__host__ __device__ __inline__ void invRotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real distX = posX - originX;
    real distY = posY - originY;
    real distZ = posZ - originZ;

    newPosX = distX;
    newPosY =  distY*cos(angle) + distZ*sin(angle);
    newPosZ = -distY*sin(angle) + distZ*cos(angle);

    newPosX += originX;
    newPosY += originY;
    newPosZ += originZ;
}

__host__ __device__ __inline__ void rotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real distX = posX - originX;
    real distY = posY - originY;
    real distZ = posZ - originZ;

    newPosX =  distX*cos(angle) + distZ*sin(angle);
    newPosY =  distY;
    newPosZ = -distX*sin(angle) + distZ*cos(angle);

    newPosX += originX;
    newPosY += originY;
    newPosZ += originZ;
}

__host__ __device__ __inline__ void invRotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real distX = posX - originX;
    real distY = posY - originY;
    real distZ = posZ - originZ;

    newPosX = distX*cos(angle) - distZ*sin(angle);
    newPosY = distY;
    newPosZ = distX*sin(angle) + distZ*cos(angle);

    newPosX += originX;
    newPosY += originY;
    newPosZ += originZ;
}

__host__ __device__ __inline__ void rotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real distX = posX - originX;
    real distY = posY - originY;
    real distZ = posZ - originZ;

    newPosX = distX*cos(angle) - distY*sin(angle);
    newPosY = distX*sin(angle) + distY*cos(angle);
    newPosZ = distZ;

    newPosX += originX;
    newPosY += originY;
    newPosZ += originZ;
}

__host__ __device__ __inline__ void invRotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real distX = posX - originX;
    real distY = posY - originY;
    real distZ = posZ - originZ;

    newPosX =  distX*cos(angle) + distY*sin(angle);
    newPosY = -distX*sin(angle) + distY*cos(angle);
    newPosZ = distZ;

    newPosX += originX;
    newPosY += originY;
    newPosZ += originZ;
}

__host__ __device__ __inline__ void rotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY)
{
    newPosX = posX*cos(angle) - posY*sin(angle);
    newPosY = posX*sin(angle) + posY*cos(angle);  
}

__host__ __device__ __inline__ void invRotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY)
{
    newPosX =  posX*cos(angle) + posY*sin(angle);
    newPosY = -posX*sin(angle) + posY*cos(angle);  
}

__host__ __device__ __inline__ void rotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX = posX;
    newPosY = posY*cos(angle) - posZ*sin(angle);
    newPosZ = posY*sin(angle) + posZ*cos(angle);
}

__host__ __device__ __inline__ void invRotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX = posX;
    newPosY =  posY*cos(angle) + posZ*sin(angle);
    newPosZ = -posY*sin(angle) + posZ*cos(angle);
}

__host__ __device__ __inline__ void rotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX =  posX*cos(angle) + posZ*sin(angle);
    newPosY =  posY;
    newPosZ = -posX*sin(angle) + posZ*cos(angle);
}

__host__ __device__ __inline__ void invRotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX = posX*cos(angle) - posZ*sin(angle);
    newPosY = posY;
    newPosZ = posX*sin(angle) + posZ*cos(angle);
}

__host__ __device__ __inline__ void rotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX = posX*cos(angle) - posY*sin(angle);
    newPosY = posX*sin(angle) + posY*cos(angle);
    newPosZ = posZ;
}

__host__ __device__ __inline__ void invRotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX =  posX*cos(angle) + posY*sin(angle);
    newPosY = -posX*sin(angle) + posY*cos(angle);
    newPosZ = posZ;
}

__host__ __device__ uint find_nearest_cellBSW(uint index, 
                                              real* coordsX, real* coordsY, real* coordsZ, 
                                              real posX, real posY, real posZ, 
                                              uint* neighborsX, uint* neighborsY, uint* neighborsZ, uint* neighborsWSB){   
    uint new_index = index;

    while(coordsX[new_index] > posX && coordsY[new_index] > posY && coordsZ[new_index] > posZ ){ new_index = max(1, neighborsWSB[new_index]);}

    while(coordsX[new_index] > posX && coordsY[new_index] > posY ){ new_index = max(1, neighborsZ[neighborsWSB[new_index]]);}
    while(coordsX[new_index] > posX && coordsZ[new_index] > posZ ){ new_index = max(1, neighborsY[neighborsWSB[new_index]]);}
    while(coordsY[new_index] > posY && coordsZ[new_index] > posZ ){ new_index = max(1, neighborsX[neighborsWSB[new_index]]);}

    while(coordsX[new_index] > posX){ new_index = max(1, neighborsY[neighborsZ[neighborsWSB[new_index]]]);}
    while(coordsY[new_index] > posY){ new_index = max(1, neighborsX[neighborsZ[neighborsWSB[new_index]]]);}
    while(coordsZ[new_index] > posZ){ new_index = max(1, neighborsX[neighborsY[neighborsWSB[new_index]]]);}

    while(coordsX[new_index] < posX){ new_index = max(1, neighborsX[new_index]);}
    while(coordsY[new_index] < posY){ new_index = max(1, neighborsY[new_index]);}
    while(coordsZ[new_index] < posZ){ new_index = max(1, neighborsZ[new_index]);}

    return neighborsWSB[new_index];
}






#endif
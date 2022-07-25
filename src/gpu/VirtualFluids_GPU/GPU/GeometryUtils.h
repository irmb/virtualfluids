#ifndef _GEOMETRYUTILS_H
#define _GEOMETRYUTILS_H

__inline__ __host__ __device__ void getNeighborIndicesOfBSW(  uint k, //index of DIR_MMM node
                                        uint &ke, uint &kn, uint &kt, uint &kne, uint &kte,uint &ktn, uint &ktne,
                                        uint* neighborX, uint* neighborY, uint* neighborZ)
{
    ke   = neighborX[k];
    kn   = neighborY[k];
    kt   = neighborZ[k];
    kne  = neighborY[ke];
    kte  = neighborZ[ke];
    ktn  = neighborZ[kn];
    ktne = neighborX[ktn];
}

__inline__ __host__ __device__ uint findNearestCellBSW(uint index, 
                                              real* coordsX, real* coordsY, real* coordsZ, 
                                              real posX, real posY, real posZ, 
                                              uint* neighborsX, uint* neighborsY, uint* neighborsZ, uint* neighborsWSB)
{
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

__inline__ __host__ __device__ void getInterpolationWeights(real &dW, real &dE, real &dN, real &dS, real &dT, real &dB,
                                        real tmpX, real tmpY, real tmpZ)
{
    dW = tmpX;      
    dE = 1.f - dW;        
    dS = tmpY;    
    dN = 1.f - dS;      
    dB = tmpZ;         
    dT = 1.f - dB;     
}

__inline__ __host__ __device__ real trilinearInterpolation( real dW, real dE, real dN, real dS, real dT, real dB,
                                        uint k,  uint ke, uint kn, uint kt, uint kne, uint kte, uint ktn, uint ktne,
                                        real* quantity )
{
    real interpolatedValue = (            dE*dN*dT*quantity[k]    + dW*dN*dT*quantity[ke]
                                        + dE*dS*dT*quantity[kn]   + dW*dS*dT*quantity[kne]
                                        + dE*dN*dB*quantity[kt]   + dW*dN*dB*quantity[kte]
                                        + dE*dS*dB*quantity[ktn]  + dW*dS*dB*quantity[ktne] );
    return interpolatedValue;
}

__inline__ __host__ __device__ void translate2D(real &posX, real &posY, real &newPosX, real &newPosY, real &translationX, real &translationY)
{
    newPosX = posX + translationX;
    newPosY = posY + translationY;
}

__inline__ __host__ __device__ void invTranslate2D(real &posX, real &posY, real &newPosX, real &newPosY, real &translationX, real &translationY)
{
    newPosX = posX - translationX;
    newPosY = posY - translationY;
}

__inline__ __host__ __device__ void translate3D(real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &translationX, real &translationY, real &translationZ)
{
    newPosX = posX + translationX;
    newPosY = posY + translationY;
    newPosZ = posZ + translationZ;
}

__inline__ __host__ __device__ void invTranslate3D(real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &translationX, real &translationY, real &translationZ)
{
    newPosX = posX - translationX;
    newPosY = posY - translationY;
    newPosZ = posZ - translationZ;
}

__inline__ __host__ __device__ void rotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY)
{
    newPosX = posX*cos(angle) - posY*sin(angle);
    newPosY = posX*sin(angle) + posY*cos(angle);  
}

__inline__ __host__ __device__ void rotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY, real &originX, real &originY)
{
    real tmpX, tmpY;
    invTranslate2D(posX, posY, newPosX, newPosY, originX, originY);
    rotate2D(angle, newPosX, newPosY, tmpX, tmpY); 
    translate2D(tmpX, tmpY, newPosX, newPosY, originX, originY);
}

__inline__ __host__ __device__ void invRotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY)
{
    newPosX =  posX*cos(angle) + posY*sin(angle);
    newPosY = -posX*sin(angle) + posY*cos(angle);  
}

__inline__ __host__ __device__ void invRotate2D(real &angle, real &posX, real &posY, real &newPosX, real &newPosY, real &originX, real &originY)
{
    real tmpX, tmpY;
    invTranslate2D(posX, posY, newPosX, newPosY, originX, originY);
    invRotate2D(angle, newPosX, newPosY, tmpX, tmpY); 
    translate2D(tmpX, tmpY, newPosX, newPosY, originX, originY);
}

__inline__ __host__ __device__ void rotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX = posX;
    rotate2D(angle, posY, posZ, newPosY, newPosZ);
}

__inline__ __host__ __device__ void rotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    rotateAboutX3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void invRotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosX = posX;
    invRotate2D(angle, posY, posZ, newPosY, newPosZ);
}

__inline__ __host__ __device__ void invRotateAboutX3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    invRotateAboutX3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void rotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{    
    newPosY =  posY;
    rotate2D(angle, posX, posZ, newPosX, newPosZ);
}

__inline__ __host__ __device__ void rotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    rotateAboutY3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void invRotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosY =  posY;
    invRotate2D(angle, posX, posZ, newPosX, newPosZ);
}

__inline__ __host__ __device__ void invRotateAboutY3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    invRotateAboutY3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}


__inline__ __host__ __device__ void rotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosZ = posZ;
    rotate2D(angle, posX, posY, newPosX, newPosY);
}

__inline__ __host__ __device__ void rotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    rotateAboutZ3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void invRotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ)
{
    newPosZ = posZ;
    invRotate2D(angle, posX, posY, newPosX, newPosY);
}

__inline__ __host__ __device__ void invRotateAboutZ3D(real &angle, real &posX, real &posY, real &posZ, real &newPosX, real &newPosY, real &newPosZ, real &originX, real &originY, real &originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    invRotateAboutZ3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}


#endif
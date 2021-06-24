#ifndef _GEOMETRYUTILS_H
#define _GEOMETRYUTILS_H

__inline__ __host__ __device__ void getNeighborIndicesBSW(  uint k, //index of BSW node
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

__inline__ __host__ __device__ void getInterpolationWeights(real &dW, real &dE, real &dN, real &dS, real &dT, real &dB,
                                        real distX, real distY, real distZ)
{
    dW = distX;      
    dE = 1.f - dW;        
    dS = distY;    
    dN = 1.f - dS;      
    dB = distZ;         
    dT = 1.f - dB;     
}

__inline__ __host__ __device__ real trilinearInterpolation( real dW, real dE, real dN, real dS, real dT, real dB,
                                        uint k,  uint ke, uint kn, uint kt, uint kne, uint kte,uint ktn, uint ktne,
                                        real* quantity )
{
    real interpolatedValue = 0.125f*(     dE*dN*dT*quantity[k]    + dW*dN*dT*quantity[ke]
                                        + dE*dS*dT*quantity[kn]   + dW*dS*dT*quantity[kne]
                                        + dE*dN*dB*quantity[kt]   + dW*dN*dB*quantity[kte]
                                        + dE*dS*dB*quantity[ktn]  + dW*dS*dB*quantity[ktne] );
    return interpolatedValue;
}


#endif
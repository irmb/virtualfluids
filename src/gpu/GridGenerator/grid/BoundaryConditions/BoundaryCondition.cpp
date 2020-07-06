#include "BoundaryCondition.h"

#include <cmath>

#include "grid/BoundaryConditions/Side.h"
#include "grid/Grid.h"

bool BoundaryCondition::isSide( SideType side ) const
{
    return this->side->whoAmI() == side;
}

VF_PUBLIC void VelocityBoundaryCondition::setVelocityProfile(SPtr<Grid> grid, std::function<void(real, real, real, real&, real&, real&)> velocityProfile)
{
    for( uint index = 0; index < this->indices.size(); index++ ){

            real x, y, z;

            grid->transIndexToCoords( this->indices[index], x, y, z );

            velocityProfile(x,y,z,this->vxList[index],this->vyList[index],this->vzList[index]);
    }
}

void GeometryBoundaryCondition::setTangentialVelocityForPatch(SPtr<Grid> grid, uint patch, 
                                                              real p1x, real p1y, real p1z, 
                                                              real p2x, real p2y, real p2z, 
                                                              real v, real r)
{
    for( uint index = 0; index < this->indices.size(); index++ ){
        if( this->patches[index] == patch ){

            real x, y, z;

            grid->transIndexToCoords( this->indices[index], x, y, z );

            real pVecX = p2x - p1x;
            real pVecY = p2y - p1y;
            real pVecZ = p2z - p1z;

            real xVecX =   x - p1x;
            real xVecY =   y - p1y;
            real xVecZ =   z - p1z;

            // compute and normalize tangent

            real tangentX = pVecY * xVecZ - pVecZ * xVecY;
            real tangentY = pVecZ * xVecX - pVecX * xVecZ; 
            real tangentZ = pVecX * xVecY - pVecY * xVecX;

            real tangentNorm = sqrt( tangentX*tangentX + tangentY*tangentY + tangentZ*tangentZ );

            tangentX /= tangentNorm;
            tangentY /= tangentNorm;
            tangentZ /= tangentNorm;

            // compute distance from rotation axis

            real projection = ( pVecX*xVecX + pVecY*xVecY + pVecZ*xVecZ ) / ( pVecX*pVecX + pVecY*pVecY + pVecZ*pVecZ );

            real d = sqrt( ( xVecX - projection * pVecX ) * ( xVecX - projection * pVecX )
                         + ( xVecY - projection * pVecY ) * ( xVecY - projection * pVecY )
                         + ( xVecZ - projection * pVecZ ) * ( xVecZ - projection * pVecZ ) );

            this->vxList[index] = tangentX * d / r * v;
            this->vyList[index] = tangentY * d / r * v;
            this->vzList[index] = tangentZ * d / r * v;

        }
    }
}

#ifndef DataDefinitions_H
#define DataDefinitions_H

#include "PassiveScalar.h"

#define SOA

//////////////////////////////////////////////////////////////////////////

#define LENGTH_VECTOR       3

#ifdef USE_PASSIVE_SCALAR
    #define LENGTH_CELL_DATA    6
#else
    #define LENGTH_CELL_DATA    5
#endif

#define LENGTH_CELL_TO_CELL 6

#define LENGTH_FACE_TO_CELL 2

#define LENGTH_FINE_TO_COARSE 9

#define LENGTH_COARSE_TO_FINE 9

//////////////////////////////////////////////////////////////////////////

#ifdef SOA

#define VEC_X(vecIdx, numberOfVectors)  ( 0 * numberOfVectors + vecIdx )
#define VEC_Y(vecIdx, numberOfVectors)  ( 1 * numberOfVectors + vecIdx )
#define VEC_Z(vecIdx, numberOfVectors)  ( 2 * numberOfVectors + vecIdx )
                                                           
#define RHO__( cellIdx, numberOfCells ) ( 0 * numberOfCells   + cellIdx )
#define RHO_U( cellIdx, numberOfCells ) ( 1 * numberOfCells   + cellIdx )
#define RHO_V( cellIdx, numberOfCells ) ( 2 * numberOfCells   + cellIdx )
#define RHO_W( cellIdx, numberOfCells ) ( 3 * numberOfCells   + cellIdx )
#define RHO_E( cellIdx, numberOfCells ) ( 4 * numberOfCells   + cellIdx )

#ifdef USE_PASSIVE_SCALAR
    #define RHO_S( cellIdx, numberOfCells ) ( 5 * numberOfCells   + cellIdx )
#endif // USE_PASSIVE_SCALAR

#define CELL_TO_CELL( cellIdx, neighborIdx, numberOfCells ) ( neighborIdx * numberOfCells + cellIdx )

#define NEG_CELL( faceIdx, numberOfFaces ) (                 faceIdx )
#define POS_CELL( faceIdx, numberOfFaces ) ( numberOfFaces + faceIdx )

#define FINE_TO_COARSE( idx, cellIdx, number ) ( cellIdx * number + idx )
#define COARSE_TO_FINE( idx, cellIdx, number ) ( cellIdx * number + idx )

#endif

//////////////////////////////////////////////////////////////////////////

#ifdef AOS

#define VEC_X(vecIdx, numberOfVectors)  ( vecIdx * LENGTH_VECTOR     )
#define VEC_Y(vecIdx, numberOfVectors)  ( vecIdx * LENGTH_VECTOR + 1 )
#define VEC_Z(vecIdx, numberOfVectors)  ( vecIdx * LENGTH_VECTOR + 2 )
                                                           
#define RHO__( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA     )
#define RHO_U( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 1 )
#define RHO_V( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 2 )
#define RHO_W( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 3 )
#define RHO_E( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 4 )

#ifdef USE_PASSIVE_SCALAR
    #define RHO_S( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 5 )
#endif // USE_PASSIVE_SCALAR
                                                                         
#define CELL_TO_CELL( cellIdx, neighborIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_TO_CELL + neighborIdx )

#define NEG_CELL( faceIdx, numberOfFaces ) ( faceIdx * LENGTH_FACE_TO_CELL     )
#define POS_CELL( faceIdx, numberOfFaces ) ( faceIdx * LENGTH_FACE_TO_CELL + 1 )

#define FINE_TO_COARSE( idx, cellIdx, number ) ( cellIdx * LENGTH_FINE_TO_COARSE + idx )
#define COARSE_TO_FINE( idx, cellIdx, number ) ( cellIdx * LENGTH_COARSE_TO_FINE + idx )

#endif

//////////////////////////////////////////////////////////////////////////

#endif

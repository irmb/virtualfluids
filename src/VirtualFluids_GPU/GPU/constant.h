#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#include <cuda_runtime.h>
#include "LBM/LB.h"

#ifdef ISDOUBLE
#define c1o2		0.5
#define c3o2		1.5
#define c1o3		0.333333333333333
#define c2o3		0.666666666666667
#define c1o4		0.25
#define c3o4		0.75
#define c1o6		0.166666666666667
#define c1o7		0.142857142857143
#define c1o8		0.125
#define c1o9		0.111111111111111
#define c2o9		0.222222222222222
#define c4o9		0.444444444444444
#define c1o10		0.1
#define c1o12		0.083333333333333
#define c1o16		0.0625
#define c3o16		0.1875
#define c9o16		0.5625
#define c1o18		0.055555555555556
#define c1o20		0.05
#define c19o20		0.95
#define c21o20		1.05
#define c1o24		0.041666666666667
#define c1o27		0.037037037037037
#define c3o32		0.09375
#define c4o32		0.125
#define c1o36		0.027777777777778
#define c1o48		0.020833333333333
#define c1o64		0.015625
#define c3o64		0.046875
#define c9o64		0.140625
#define c27o64		0.421875
#define c1o66		0.015151515151515
#define c1o72		0.013888888888889
#define c1o264		0.003787878787879
#define c8over27	0.296296296296296
#define c2over27	0.074074074074074
#define c1over54	0.018518518518519
#define c1o100		0.01
#define c99o100		0.99
#define c1over126	0.007936507936508
#define c1over216	0.004629629629630
#define c9over2		4.5

#define zero		0.
#define one			1.
#define two			2.
#define three		3.
#define four		4.
#define five		5.
#define six			6.
#define seven		7.
#define eight		8.
#define nine		9.
#define ten 		10.
#define eleven  	11.
#define twelve  	12.
#define fourteen	14.
#define fiveteen	15.
#define sixteen 	16.
#define seventeen 	17.
#define eighteen	18.
#define twentyone	21.
#define twentyfour	24.
#define twentyseven	27.
#define twentyeight	28.
#define twentynine	29.
#define thirty  	30.
#define thirtytwo	32.
#define thirtythree	33.
#define thirtyfour	34.
#define thirtysix	36.
#define fourty		40.
#define fortytwo	42.
#define fourtyeight	48.
#define fifty		50.
#define fiftyfour	54.
#define fiftysix	56.
#define sixtyfour	64.
#define sixtysix	66.
#define sixtyeight	68.
#define seventytwo	72.
#define eightyfour	84.
#define ninetysix	96.
#define c130		130.0
#define c152		152.0
#define c166		166.0
#define c195		195.0
#define c290		290.0
#define c367		367.0

#define Op0000002   0.0000002
#define smallSingle 0.0000000002f
#endif

#ifdef ISFLOAT
#define c1o2		0.5f
#define c3o2		1.5f
#define c1o3		0.333333f
#define c2o3		0.666667f
#define c1o4		0.25f
#define c3o4		0.75f
#define c1o6		0.166667f
#define c1o7		0.142857f
#define c1o8		0.125f
#define c1o9		0.111111f
#define c2o9		0.222222f
#define c4o9		0.444444f
#define c1o10		0.1f
#define c1o12		0.083333f
#define c1o16		0.0625f
#define c3o16		0.1875f
#define c9o16		0.5625f
#define c1o18		0.055556f
#define c1o20		0.05f
#define c19o20		0.95f
#define c21o20		1.05f
#define c1o24		0.041667f
#define c1o27		0.037037f
#define c3o32		0.09375f
#define c4o32		0.125f
#define c1o36		0.027778f
#define c1o48		0.020833f
#define c1o64		0.015625f
#define c3o64		0.046875f
#define c9o64		0.140625f
#define c27o64		0.421875f
#define c1o66		0.015152f
#define c1o72		0.013889f
#define c1o264		0.003788f
#define c8over27	0.296296f
#define c2over27	0.074074f
#define c1over54	0.018519f
#define c1o100		0.01f
#define c99o100		0.99f
#define c1over126	0.007937f
#define c1over216	0.004630f
#define c9over2		4.5f

#define zero		0.f
#define one			1.f
#define two			2.f
#define three		3.f
#define four		4.f
#define five		5.f
#define six			6.f
#define seven		7.f
#define eight		8.f
#define nine		9.f
#define ten 		10.f
#define eleven  	11.f
#define twelve  	12.f
#define fourteen	14.f
#define fiveteen	15.f
#define sixteen 	16.f
#define seventeen 	17.f
#define eighteen	18.f
#define twentyone	21.f
#define twentyfour	24.f
#define twentyseven	27.f
#define twentyeight	28.f
#define twentynine	29.f
#define thirty  	30.f
#define thirtytwo	32.f
#define thirtythree	33.f
#define thirtyfour	34.f
#define thirtysix	36.f
#define fourty		40.f
#define fortytwo	42.f
#define fourtyeight	48.f
#define fifty		50.f
#define fiftyfour	54.f
#define fiftysix	56.f
#define sixtyfour	64.f
#define sixtysix	66.f
#define sixtyeight	68.f
#define seventytwo	72.f
#define eightyfour	84.f
#define ninetysix	96.f
#define c130		130.0f
#define c152		152.0f
#define c166		166.0f
#define c195		195.0f
#define c290		290.0f
#define c367		367.0f

#define Op0000002   0.0000002f
#define smallSingle 0.0000000002f
#endif






//__constant__ doubflo c1o2			= 0.5				 ;
//__constant__ doubflo c3o2			= 1.5				 ;
//__constant__ doubflo c1o3			= 0.333333333333333	 ;
//__constant__ doubflo c2o3			= 0.666666666666667	 ;
//__constant__ doubflo c1o4			= 0.25				 ;
//__constant__ doubflo c1o6			= 0.166666666666667	 ;
//__constant__ doubflo c1o8			= 0.125				 ;
//__constant__ doubflo c1o9			= 0.111111111111111	 ;
//__constant__ doubflo c2o9			= 0.222222222222222	 ;
//__constant__ doubflo c4o9			= 0.444444444444444	 ;
//__constant__ doubflo c1o10			= 0.1				 ;
//__constant__ doubflo c1o12			= 0.083333333333333	 ;
//__constant__ doubflo c1o16			= 0.0625			 ;
//__constant__ doubflo c3o16			= 0.1875			 ;
//__constant__ doubflo c9o16			= 0.5625			 ;
//__constant__ doubflo c1o18			= 0.055555555555556	 ;
//__constant__ doubflo c1o20			= 0.05				 ;
//__constant__ doubflo c19o20			= 0.95				 ;
//__constant__ doubflo c21o20			= 1.05				 ;
//__constant__ doubflo c1o24			= 0.041666666666667	 ;
//__constant__ doubflo c1o27			= 0.037037037037037	 ;
//__constant__ doubflo c4o32			= 0.125				 ;
//__constant__ doubflo c1o36			= 0.027777777777778	 ;
//__constant__ doubflo c1o48			= 0.020833333333333	 ;
//__constant__ doubflo c1o64			= 0.015625			 ;
//__constant__ doubflo c3o64			= 0.046875			 ;
//__constant__ doubflo c9o64			= 0.140625			 ;
//__constant__ doubflo c27o64			= 0.421875			 ;
//__constant__ doubflo c1o72			= 0.013888888888889	 ;
//__constant__ doubflo c8over27		= 0.296296296296296	 ;
//__constant__ doubflo c2over27		= 0.074074074074074	 ;
//__constant__ doubflo c1over54		= 0.018518518518519	 ;
//__constant__ doubflo c1o100			= 0.01				 ;
//__constant__ doubflo c99o100		= 0.99				 ;
//__constant__ doubflo c1over126		= 0.007936507936508	 ;
//__constant__ doubflo c1over216		= 0.004629629629630	 ;
//__constant__ doubflo c9over2		= 4.5				 ;
//
//__constant__ doubflo zero			= 0.				 ;
//__constant__ doubflo one			= 1.				 ;
//__constant__ doubflo two			= 2.				 ;
//__constant__ doubflo three			= 3.				 ;
//__constant__ doubflo four			= 4.				 ;
//__constant__ doubflo five			= 5.				 ;
//__constant__ doubflo six			= 6.				 ;
//__constant__ doubflo eight			= 8.				 ;
//__constant__ doubflo nine			= 9.				 ;
//__constant__ doubflo ten 			= 10.				 ;
//__constant__ doubflo eleven  		= 11.				 ;
//__constant__ doubflo twelve  		= 12.				 ;
//__constant__ doubflo fourteen		= 14.				 ;
//__constant__ doubflo sixteen 		= 16.				 ;
//__constant__ doubflo seventeen 		= 17.				 ;
//__constant__ doubflo eighteen		= 18.				 ;
//__constant__ doubflo twentyseven	= 27.				 ;
//__constant__ doubflo thirty  		= 30.				 ;
//__constant__ doubflo thirtytwo		= 32.				 ;
//__constant__ doubflo thirtythree	= 33.				 ;
//__constant__ doubflo thirtyfour		= 34.				 ;
//__constant__ doubflo thirtysix		= 36.				 ;
//__constant__ doubflo fortytwo		= 42.				 ;
//__constant__ doubflo fiftyfour		= 54.				 ;
//__constant__ doubflo sixtyfour		= 64.				 ;
//__constant__ doubflo sixtysix		= 66.				 ;
//__constant__ doubflo sixtyeight		= 68.				 ;
//__constant__ doubflo seventytwo		= 72.				 ;
//__constant__ doubflo eightyfour		= 84.				 ;
//__constant__ doubflo ninetysix		= 96.				 ;
//
//__constant__ doubflo Op0000002		= 0.0000002			 ;

#endif

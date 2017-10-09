// TODO
// define
// c++

#ifndef __UTIL__
#define __UTIL__


#ifndef DRAWNOW
#define DRAWNOW mexEvalString("drawnow;");
#endif

#ifndef PRINTF
#define PRINTF mexPrintf
#endif

#ifndef MATLAB_COMMAND
#define MATLAB_COMMAND(cmd) mexEvalString(cmd)
#endif

#ifndef CLOGIC
#define CLOGIC 0
#endif

#ifndef MATLABLOGIC
#define MATLABLOGIC 1
#endif

#ifndef DIRECTION_X
#define DIRECTION_X 0
#endif

#ifndef DIRECTION_Y
#define DIRECTION_Y 1
#endif

#ifndef DIRECTION_Z
#define DIRECTION_Z 2
#endif



#ifndef ADD2
#define ADD2(x, y, height) ((height)*(x) + (y))
#endif

#ifndef ADD3
#define ADD3(x, y, z, width, height) ((width)*(height)*(z) + height*(x) + (y))
#endif

#ifndef ADD2d
#define ADD2d(x, y, width, height) ( (height)*(x)+(y) )
#endif

//#ifndef ADD3d
//#define ADD3d(x, y, z, width, height, depth) ( (width)*(height)*(z) + (height)*(x) + (y) )
//#endif

#ifndef ADD3d
#define ADD3d(x, y, z, width, height, depth) ( (width)*(height)*(z) + (width)*(y) + (x) )
#endif

#ifndef ADD2dC
#define ADD2dC(x, y, width, height) (width)*(y)+(x)
#endif

#ifndef ADD3dC
#define ADD3dC(x, y, z, width, height, depth) ( (width)*((height)*(z)+(y))+(x) )
#endif

#ifndef MIN
	#define MIN(x,y)		( ( x >= y ) ? y : x )
#endif

#ifndef MAX
	#define MAX(x,y)		( ( x >= y ) ? x : y )
#endif

#ifndef MIN_MAX_THRESHOLD
#define MIN_MAX_THRESHOLD(a, mn, mx) MIN(MAX(a, mn), mx)
#endif

#ifndef SGN
	#define SGN(x)			( ( x >= 0 ) ? (1) : (-1) )
#endif

// la meme chose que SGN mais en flotant
#ifndef SGNF
#define SGNF(x)			( ( x >= 0. ) ? (1) : (-1) )
#endif
    
#if !defined (PI)
	#define PI				 3.14159265358979
#endif

#if !defined (PI_DOUBLE)
	#define PI_DOUBLE  6.28318530717959
#endif

#ifndef NULL
#define NULL    ((void *)0)
#endif

#define INV3  .333333333333333f
#define INV6  .166666666666667f
#define ROOT3 1.732050807568877f
#if !defined (PI_3)
#define PI_3 1.047197551196598f
#endif


#ifndef SQR
#define SQR(x) (x)*(x)
#endif

#ifndef SWAP_TEMP
#define SWAP_TEMP(a, b, temp) {(temp) = (a); (a) = (b); (b) = (temp);}
#endif


#ifndef DIST2
#define DIST2(x1, y1, x2, y2) (sqrt(SQR(x2-x1)+SQR(y2-y1)))
#endif

#ifndef DOT2
#define DOT2(a1, a2, b1, b2) ( (a1) * (b1) + (a2) * (b2) )
#endif


#ifndef DET22
#define DET22(a11, a12, a21, a22) ( (a11)*(a22) - (a21)*(a12) )
#endif

#define MAT2_INV_A11(det, a11, a12, a21, a22) (a22 / (det))
#define MAT2_INV_A12(det, a11, a12, a21, a22) (-a12 / (det))
#define MAT2_INV_A21(det, a11, a12, a21, a22) (-a21 / (det))
#define MAT2_INV_A22(det, a11, a12, a21, a22) (a11 / (det))

#define MAT2_INV(a11, a12, a21, a22, r11, r12, r21, r22) { \
double det = DET22(a11, a12, a21, a22); \
r11 = MAT2_INV_A11(det, a11, a12, a21, a22); \
r12 = MAT2_INV_A12(det, a11, a12, a21, a22); \
r21 = MAT2_INV_A21(det, a11, a12, a21, a22); \
r22 = MAT2_INV_A22(det, a11, a12, a21, a22); }

#define MAT2_MULT_VECT(a11, a12, a21, a22, v1, v2, r1, r2) { \
(r1) = DOT2(a11, a12, v1, v2); \
(r2) = DOT2(a21, a22, v1, v2); }

//#ifndef det33
//#define det33(a11, a12, a13, a21, a22, a23, a31, a32, a33) (((a11)*(a22)*(a33)+(a21)*(a32)*(a13)+(a31)*(a12)*(a23)) - ((a31)*(a22)*(a13)+(a21)*(a12)*(a33)+(a11)*(a32)*(a23)))
//#endif

#ifndef DET33
#define DET33(a11, a12, a13, a21, a22, a23, a31, a32, a33) ( (a11) * ( (a22) * (a33) - (a32) * (a23) ) - (a12) * ( (a21) * (a33) - (a23) * (a31) ) + (a13) * ( (a21) * (a32) - (a22) * (a31)) )
#endif

#define DET33_V2(A) DET33(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8])

#ifndef DOT3
#define DOT3(a1, a2, a3, b1, b2, b3) ( (a1) * (b1) + (a2) * (b2) + (a3) * (b3) )
#endif

#define TRANSPOSE33(a11, a12, a13, a21, a22, a23, a31, a32, a33) { \
 double temp; \
 SWAP_TEMP(a12, a21, temp)\
 SWAP_TEMP(a13, a31, temp)\
 SWAP_TEMP(a23, a32, temp) }

#define TRANSPOSE33_V2(A) TRANSPOSE33(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8])

#define MAT3_INV_A11(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) (  DET22(a22, a23, a32, a33) / (det) )
#define MAT3_INV_A12(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) ( -DET22(a12, a13, a32, a33) / (det) )
#define MAT3_INV_A13(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) (  DET22(a12, a13, a22, a23) / (det) )
#define MAT3_INV_A21(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) ( -DET22(a21, a23, a31, a33) / (det) )
#define MAT3_INV_A22(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) (  DET22(a11, a13, a31, a33) / (det) )
#define MAT3_INV_A23(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) ( -DET22(a11, a13, a21, a23) / (det) )
#define MAT3_INV_A31(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) (  DET22(a21, a22, a31, a32) / (det) )
#define MAT3_INV_A32(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) ( -DET22(a11, a12, a31, a32) / (det) )
#define MAT3_INV_A33(det, a11, a12, a13, a21, a22, a23, a31, a32, a33) (  DET22(a11, a12, a21, a22) / (det) )

#define MAT3_INV(a11, a12, a13, a21, a22, a23, a31, a32, a33, r11, r12, r13, r21, r22, r23, r31, r32, r33) { \
	double det = DET33(a11, a12, a13, a21, a22, a23, a31, a32, a33); \
	/*if ( det <= 0.0 ) det = 0.001; \*/ \
	/* det = 1000.0; */ \
	r11 = MAT3_INV_A11(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r12 = MAT3_INV_A12(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r13 = MAT3_INV_A13(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r21 = MAT3_INV_A21(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r22 = MAT3_INV_A22(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r23 = MAT3_INV_A23(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r31 = MAT3_INV_A31(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r32 = MAT3_INV_A32(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  \
	r33 = MAT3_INV_A33(det, a11, a12, a13, a21, a22, a23, a31, a32, a33);  }

#define MAT3_INV_V2(A, R) MAT3_INV(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8], R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8])

#define MAT3_MULT_VECT(a11, a12, a13, a21, a22, a23, a31, a32, a33, v1, v2, v3, r1, r2, r3) {\
(r1) = DOT3(a11, a12, a13, v1, v2, v3); \
(r2) = DOT3(a21, a22, a23, v1, v2, v3); \
(r3) = DOT3(a31, a32, a33, v1, v2, v3); }

#define MAT3_MULT_VECT_V2(A, V, R) MAT3_MULT_VECT(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8], V[0], V[1], V[2], R[0], R[1], R[2])

#define MAT3_MULT_VECT_V3(A, v1, v2, v3, r1, r2, r3) MAT3_MULT_VECT(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8], v1, v2, v3, r1, r2, r3)



#define BILINEAR_INTERP(v00, v01, v10, v11, x, y) ((1.0-(y)) * (v00*(1.0-(x))+v01*(x)) + (y) * (v11*(x)+v10*(1.0-(x))))
#define LINEAR_INTERP(v00, v11, x) ( ((v11)-(v00))*(x) + v00 )

#ifndef NORM2
#define NORM2(x1, x2) sqrt(SQR(x1)+SQR(x2))
#endif

#ifndef NORM3
#define NORM3(x1, x2, x3) sqrt(SQR(x1)+SQR(x2)+SQR(x3))
#endif


#ifndef CALLOC
#define CALLOC(size, format) (format*)calloc(size, sizeof(format))
#endif


#ifndef FREE
#define FREE(x){ \
if ( x != NULL ) { \
        free(x); \
        x = NULL;} }
#endif

#ifndef CALLOC_TAB
#define CALLOC_TAB(ptr, nb, size, element)	{ \
												long i; \
												ptr = (element **)calloc(nb, sizeof(element *)); \
												for (i=0; i<nb; i++) (ptr)[i] = (element *)calloc(size, sizeof(element)); \
											}
#endif

#ifndef FREE_TAB
#define FREE_TAB(ptr, nb)   { \
								long i; \
								if ( ptr != NULL ) for (i=0; i<nb; i++) FREE((ptr)[i]) \
								FREE(ptr) \
							}
#endif
    

#ifndef CALLOC_VTAB
#define CALLOC_VTAB(ptr, nb, size, element)	{ \
												long i; \
												element * tmp = NULL; \
												tmp = (element *)calloc(nb*size, sizeof(element)); \
												ptr = (element **)calloc(nb, sizeof(element *));\
												for (i=0; i<nb; i++) (ptr)[i] = tmp+i*size; \
											}
#endif

#ifndef FREE_VTAB
#define FREE_VTAB(ptr)   { \
								if ( ptr != NULL ) FREE(ptr[0]);\
								FREE(ptr) \
							}    
    
#endif


#ifndef SAWP
#define SWAP(_a_, _b_, _temp_) { \
	_temp_ = _b_; \
	_b_ = _a_; \
	_a_ = _temp_; }
#endif


#define VECTOR_NORMALIZE(v1, v2, v3) { \
	float norm; \
	norm = sqrtf(v1*v1 + v2*v2 + v3*v3); \
	if ( norm != 0.0f ) { \
		v1 /= norm; \
		v2 /= norm; \
		v3 /= norm; } \
else { v1 = v3 = 0.0f; v2 = -1.0f; } }

#define CROSS(a1, a2, a3, b1, b2, b3, r1, r2, r3){ \
	r1 = (a2)*(b3)-(a3)*(b2); \
	r2 = (a3)*(b1)-(a1)*(b3); \
	r3 = (a1)*(b2)-(a2)*(b1); }

#define SPROD_VMV(v1, v2, v3, m11, m12, m13, m21, m22, m23, m31, m32, m33, w1, w2, w3) \
	( (v1)*((m11)*(w1)+(m12)*(w2)+(m13)*(w3)) + (v2)*((m21)*(w1)+(m22)*(w2)+(m23)*(w3)) + (v3)*((m31)*(w1)+(m32)*(w2)+(m33)*(w3)) )
























#define KMAX .001f
#define INV_24 .041666666666f
#define INV_4 .25f

#define VAR_THRESHOLD(v, th) { if ((v) > (th)) (v) = (th); else if ((v) < -(th)) (v) = -(th); }

#define NORMAL_THRESHOLD(vx, vy, vz, th) { \
	float dipx, dipz; \
	if ( vy != 0.0f ) \
		{ \
		dipx = -vx / vy; \
		dipz = -vz / vy; \
		} \
	  else \
		{ \
	  dipx = 0.0f; \
	  dipz = 0.0f; \
		} \
	VAR_THRESHOLD(dipx, th) \
	VAR_THRESHOLD(dipz, th) \
	vy = -1.0f; \
	vx = -dipx; \
	vz = -dipz; \
	VECTOR_NORMALIZE(vx, vy, vz) }

#define QUADRIC_YCOEF(K1, v1x, v1y, v1z, K2, v2x, v2y, v2z, nx, ny, nz, y0, y1, y2, y3, y4, y5, y6, y7, y8) { \
	y0 = K1 / 2.0f * v1x * v1x + K2 / 2.0f * v2x * v2x; \
	y1 = K1 / 2.0f * v1y * v1y + K2 / 2.0f * v2y * v2y; \
	y2 = K1 / 2.0f * v1z * v1z + K2 / 2.0f * v2z * v2z; \
	y3 = K1 / 2.0f * v1x * v1y + K2 / 2.0f * v2x * v2y; \
	y4 = K1 / 2.0f * v1x * v1z + K2 / 2.0f * v2x * v2z; \
	y5 = K1 / 2.0f * v1y * v1z + K2 / 2.0f * v2y * v2z; \
	y6 = - nx; \
	y7 = - ny; \
	y8 = - nz; }

#define TAYLOR_CURVATURE(Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, A0, A1, A2, A3, A4, th) { \
	NORMAL_THRESHOLD(Y6, Y7, Y8, th) \
	if ( Y7 != 0.0f ) \
			{ \
		A0 = -Y6 / Y7; \
		A1 = -Y8 / Y7; \
		A2 = -(Y0+2.0f*Y3*A0+Y1*A0*A0)/Y7; \
		A3 = -2.0f*(Y4+Y5*A0+Y3*A1+Y1*A0*A1)/Y7; \
		A4 = -(Y2+2.0f*Y5*A1+Y1*A1*A1)/Y7; \
			} \
				else \
				{ \
		A0 = A1 = A2 = A3 = A4 = 0.0f; \
				} }


#define THREADS_START(N, attr, th_id, funct, data) { \
	int iithreads; \
	for (iithreads=0; iithreads<N; iithreads++) \
	{ \
		pthread_attr_init (&attr[iithreads]); \
		pthread_attr_setdetachstate (&attr[iithreads], PTHREAD_CREATE_JOINABLE); \
		pthread_create (&th_id[iithreads], &attr[iithreads], funct, &data[iithreads]); } }

#define THREADS_JOIN(N, th_id) { \
	int iithreads; \
	for (iithreads=0; iithreads<N; iithreads++) \
	pthread_join(th_id[iithreads], NULL); }

#define ARRAY_CHAR_COPY(_src, _dst) {\
  FREE(_dst) \
  if ( _src != NULL ) { _dst = (char*)calloc(strlen(_src) + 1, sizeof(char)); strcpy(_dst, _src); } \
}



#define EPS_EQUAL_DOUBLE 0.00001

#ifdef __cplusplus
extern "C" {
#endif


  long double_equal(double v1, double v2);
  int pow2(int n);
  int pow4(int n);
  int pow8(int n);

  float **malloc_array_float(int nbre, int *size);

  void **free_array_float(void **ptr, int nbre);

#ifdef __cplusplus
}
#endif

#endif
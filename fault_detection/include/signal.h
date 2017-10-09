
#ifndef __SIGNAL__
#define __SIGNAL__

#include <config.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef QUADRILINEAR_INTERP
#define QUADRILINEAR_INTERP(v00, v01, v10, v11, x, z) ((1.0-(z)) * (v00*(1.0-(x))+v01*(x)) + (z) * (v11*(x)+v10*(1.0-(x))))
#endif

	long gaussSigma2Size(double sigma);

	double gaussSize2Sigma(long size);

  EXPORT_LIB double *gaussianMask1D(double sigma, long *size);

  EXPORT_LIB double *gaussianGradMask1D(double sigma, long *size);

	EXPORT_LIB float *gaussianMask1D_float(double sigma, long *size);

	EXPORT_LIB float *gaussianGradMask1D_float(double sigma, long *size);
	
	double *gaussianMask2D(double sigmax, double sigmay, long *sizex, long *sizey);

	double *gaussianMask3D(double sigmax, double sigmay, double sigmaz, long *sizex, long *sizey, long *sizez);

	EXPORT_LIB void *gaussianMaskFree(void *mask);

	void conv1d(double *data, long size, double *mask, long sizem, double *res);
	void conv1d_float(float *data, long size, float *mask, long sizem, float *res);

	void data3d_conv_separable(double *data, long width, long height, long depth,
		double *maskx, long sizex,
		double *masky, long sizey,
		double *maskz, long sizez,
		double *res);

  EXPORT_LIB void linspace_fill(double x1, double x2, long N, double *v);


#ifdef __cplusplus
}
#endif

#endif
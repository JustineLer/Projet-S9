
#ifndef __IMG_UTILS__
#define __IMG_UTILS__

#ifdef __cplusplus
extern "C" {
#endif

#include <config.h>

	void data3d_extractx(double *data, long width, long height, long depth, long y0, long z0, double *f);
	void data3d_extracty(double *data, long width, long height, long depth, long x0, long z0, double *f);
	void data3d_extractz(double *data, long width, long height, long depth, long x0, long y0, double *f);
	void data3d_insertx(double *data, long width, long height, long depth, long y0, long z0, double *f);
	void data3d_inserty(double *data, long width, long height, long depth, long x0, long z0, double *f);
	void data3d_insertz(double *data, long width, long height, long depth, long x0, long y0, double *f);
	void data3d_extractx_float(float *data, long width, long height, long depth, long y0, long z0, float *f);
	void data3d_extracty_float(float *data, long width, long height, long depth, long x0, long z0, float *f);	
	void data3d_extractz_float(float *data, long width, long height, long depth, long x0, long y0, float *f);	
	void data3d_insertx_float(float *data, long width, long height, long depth, long y0, long z0, float *f);	
	void data3d_inserty_float(float *data, long width, long height, long depth, long x0, long z0, float *f);	
	void data3d_insertz_float(float *data, long width, long height, long depth, long x0, long y0, float *f);

#ifdef __cplusplus
}
#endif

#endif
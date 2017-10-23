
#include <stdlib.h>

#include <util.h>
#include <img_utils.h>


#define DATA3D_EXTRACTX(data, width, height, depth, y0, z0, f) {\
		long x; \
		\
    for (x = 0; x<width; x++) \
      f[x] = data[ADD3dC(x, y0, z0, width, height, depth)];\
}

#define DATA3D_EXTRACTY(data, width, height, depth, y0, z0, f) {\
	long y; \
	\
  for (y = 0; y<height; y++) \
    f[y] = data[ADD3dC(x0, y, z0, width, height, depth)]; \
}

#define DATA3D_EXTRACTZ(data, width, height, depth, y0, z0, f) {\
	long z; \
	\
  for (z = 0; z<depth; z++) \
    f[z] = data[ADD3dC(x0, y0, z, width, height, depth)]; \
}

#define DATA3D_INSERTX(data, width, height, depth, y0, z0, f) {\
	long x; \
	\
  for (x = 0; x<width; x++) \
    data[ADD3dC(x, y0, z0, width, height, depth)] = f[x]; \
}

#define DATA3D_INSERTY(data, width, height, depth, y0, z0, f) {\
	long y; \
	\
  for (y = 0; y<height; y++) \
    data[ADD3dC(x0, y, z0, width, height, depth)] = f[y]; \
}

#define DATA3D_INSERTZ(data, width, height, depth, y0, z0, f) {\
	long z; \
	\
  for (z = 0; z<depth; z++) \
    data[ADD3dC(x0, y0, z, width, height, depth)] = f[z]; \
}

void data3d_extractx(double *data, long width, long height, long depth, long y0, long z0, double *f)
{
	DATA3D_EXTRACTX(data, width, height, depth, y0, z0, f)
}

void data3d_extracty(double *data, long width, long height, long depth, long x0, long z0, double *f)
{
	DATA3D_EXTRACTY(data, width, height, depth, y0, z0, f);
}

void data3d_extractz(double *data, long width, long height, long depth, long x0, long y0, double *f)
{
	DATA3D_EXTRACTZ(data, width, height, depth, y0, z0, f)
}

void data3d_insertx(double *data, long width, long height, long depth, long y0, long z0, double *f)
{
	DATA3D_INSERTX(data, width, height, depth, y0, z0, f)
}

void data3d_inserty(double *data, long width, long height, long depth, long x0, long z0, double *f)
{
  DATA3D_INSERTY(data, width, height, depth, y0, z0, f)
}

void data3d_insertz(double *data, long width, long height, long depth, long x0, long y0, double *f)
{
	DATA3D_INSERTZ(data, width, height, depth, y0, z0, f)
}


void data3d_extractx_float(float *data, long width, long height, long depth, long y0, long z0, float *f)
{
	DATA3D_EXTRACTX(data, width, height, depth, y0, z0, f)
}

void data3d_extracty_float(float *data, long width, long height, long depth, long x0, long z0, float *f)
{
	DATA3D_EXTRACTY(data, width, height, depth, y0, z0, f);
}

void data3d_extractz_float(float *data, long width, long height, long depth, long x0, long y0, float *f)
{
	DATA3D_EXTRACTZ(data, width, height, depth, y0, z0, f)
}

void data3d_insertx_float(float *data, long width, long height, long depth, long y0, long z0, float *f)
{
	DATA3D_INSERTX(data, width, height, depth, y0, z0, f)
}

void data3d_inserty_float(float *data, long width, long height, long depth, long x0, long z0, float *f)
{
	DATA3D_INSERTY(data, width, height, depth, y0, z0, f)
}

void data3d_insertz_float(float *data, long width, long height, long depth, long x0, long y0, float *f)
{
	DATA3D_INSERTZ(data, width, height, depth, y0, z0, f)
}





//
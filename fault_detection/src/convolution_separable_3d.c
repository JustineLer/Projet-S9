
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include <util.h>
#include <signal.h>
#include <img_utils.h>
#include <convolution_separable_3d.h>

void convolution_separable_3d(double *data, long width, long height, long depth,
	double *maskx, long sizex,
	double *masky, long sizey,
	double *maskz, long sizez,
	double *res)
{
	long x, y, z, lmax;
	double *v = NULL, *vr = NULL;

	lmax = MAX(width, height);
	lmax = MAX(lmax, depth);

	v = (double*)calloc(lmax, sizeof(double));
	vr = (double*)calloc(lmax, sizeof(double));

	for (z = 0; z<depth; z++)
		for (y = 0; y<height; y++)
		{
			data3d_extractx(data, width, height, depth, y, z, v);
			conv1d(v, width, maskx, sizex, vr);
			data3d_insertx(res, width, height, depth, y, z, vr);
		}

	for (z = 0; z<depth; z++)
		for (x = 0; x<width; x++)
		{
			data3d_extracty(res, width, height, depth, x, z, v);
			conv1d(v, height, masky, sizey, vr);
			data3d_inserty(res, width, height, depth, x, z, vr);
		}

	for (y = 0; y<height; y++)
		for (x = 0; x<width; x++)
		{
			data3d_extractz(res, width, height, depth, x, y, v);
			conv1d(v, depth, maskz, sizez, vr);
			data3d_extractz(res, width, height, depth, x, y, vr);
		}
	FREE(v)
	FREE(vr)
}


// =======================================================================================
void convolutionx_3d(long type, float *in, long width, long height, long depth,
  float *mask, long mask_size,
  float *v, float *vr,
  float *out)
{
  long y, z, width2;
  float *vr2 = NULL;

  if (type == CONVOLUTION_SEPARABLE_3D_SAME)
  {
    width2 = width;
    vr2 = vr;
  }
  else if ( type == CONVOLUTION_SEPARABLE_3D_VALID)
  {
    width2 = width - mask_size + 1;
    vr2 = vr + (mask_size - 1) / 2;
  }

  for (z = 0; z<depth; z++)
    for (y = 0; y<height; y++)
    {
      data3d_extractx_float(in, width, height, depth, y, z, v);
      conv1d_float(v, width, mask, mask_size, vr);
      data3d_insertx_float(out, width2, height, depth, y, z, vr2);
    }
}

void convolutiony_3d(long type, float *in, long width, long height, long depth,
  float *mask, long mask_size,
  float *v, float *vr,
  float *out)
{
  long x, z, height2;
  float *vr2 = NULL;

  if (type == CONVOLUTION_SEPARABLE_3D_SAME)
  {
    height2 = height;
    vr2 = vr;
  }
  else if (type == CONVOLUTION_SEPARABLE_3D_VALID)
  {
    height2 = height - mask_size + 1;
    vr2 = vr + (mask_size - 1) / 2;
  }

  for (z = 0; z<depth; z++)
    for (x = 0; x<width; x++)
    {
      data3d_extracty_float(in, width, height, depth, x, z, v);
      conv1d_float(v, height, mask, mask_size, vr);
      data3d_inserty_float(out, width, height2, depth, x, z, vr2);
    }
}



void convolutionz_3d(long type, float *in, long width, long height, long depth,
  float *mask, long mask_size,
  float *v, float *vr,
  float *out)
{
  long x, y, depth2;
  float *vr2 = NULL;

  if (type == CONVOLUTION_SEPARABLE_3D_SAME)
  {
    depth2 = depth;
    vr2 = vr;
  }
  else if (type == CONVOLUTION_SEPARABLE_3D_VALID)
  {
    depth2 = depth - mask_size + 1;
    vr2 = vr + (mask_size - 1) / 2;
  }

  for (y = 0; y<height; y++)
    for (x = 0; x<width; x++)
    {
      data3d_extractz_float(in, width, height, depth, x, y, v);
      conv1d_float(v, depth, mask, mask_size, vr);
      data3d_insertz_float(out, width, height, depth2, x, y, vr2);
    }
}

void gradient_separable_valid_3d_float(long type, float *data, long width, long height, long depth,
                                 float *h_lp, float *h_hp, long h_size,
                                 float *temp1, float *temp2,
                                 float *gx, float *gy, float *gz)
{
	long lmax, width2, height2, depth2;
	float *v = NULL, *vr = NULL;

	lmax = MAX(width, height);
	lmax = MAX(lmax, depth);
  
	v = (float*)calloc(lmax, sizeof(float));
	vr = (float*)calloc(lmax, sizeof(float));

  if ( type == CONVOLUTION_SEPARABLE_3D_VALID )
  {
    width2 = width - h_size + 1;
    height2 = height - h_size + 1;
    depth2 = depth - h_size + 1;
  }
  else
  {
    width2 = width;
    height2 = height;
    depth2 = depth;
  }
  convolutionz_3d(type, data, width, height, depth, h_lp, h_size, v, vr, temp1);
  convolutiony_3d(type, temp1, width, height, depth2, h_lp, h_size, v, vr, temp2);
  convolutionx_3d(type, temp2, width, height2, depth2, h_hp, h_size, v, vr, gx);
  convolutionx_3d(type, temp1, width, height, depth2, h_lp, h_size, v, vr, temp2);
  convolutiony_3d(type, temp2, width2, height, depth2, h_hp, h_size, v, vr, gy);

  convolutionx_3d(type, data, width, height, depth, h_lp, h_size, v, vr, temp1);
  convolutiony_3d(type, temp1, width2, height, depth, h_lp, h_size, v, vr, temp2);
  convolutionz_3d(type, temp2, width2, height2, depth, h_hp, h_size, v, vr, gz);    

	FREE(v)
  FREE(vr)
}


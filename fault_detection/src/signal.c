
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include <util.h>
#include <signal.h>

long gaussSigma2Size(double sigma)
{
    return 2 * (long)ceil(3.0*sigma) + 1;
}

double gaussSize2Sigma(long size)
{
    return (double)(size-1)/(2.0*3.0);
}

double *gaussianMask1D(double sigma, long *size)
{
	double *mask = NULL, norm = 0.0, x, den;
	long width = gaussSigma2Size(sigma), i, width_2;

	width_2 = width/2;

	den = 2.0*sigma*sigma;
	mask = (double*)calloc(width, sizeof(double));
	for (i=0; i<width; i++)
	{
		x = (double)(i-width_2);
		mask[i] = exp(-(x*x)/den);
		norm += mask[i];	
	}
	for (i=0; i<width; i++)
		mask[i] /= norm;
  if ( size )
    *size = width;
	return mask;
}

double *gaussianGradMask1D(double sigma, long *size)
{
	double *mask = NULL, norm = 0.0, x, den;
	long width = gaussSigma2Size(sigma), i, width_2;

	width_2 = width/2;

	den = 2.0*sigma*sigma;
	mask = (double*)calloc(width, sizeof(double));
	for (i=0; i<width; i++)
	{
		x = (double)(i-width_2);
		mask[i] = -x*exp(-(x*x)/den);
		// norm += fabs(mask[i]);	
		norm += mask[i] * (width/2-i);
	}
	for (i=0; i<width; i++)
		mask[i] /= norm;
  if ( size )
    *size = width;
	return mask;
}


float *gaussianMask1D_float(double sigma, long *size)
{
	float *mask = NULL;
	double norm = 0.0, x, den;
	long width = gaussSigma2Size(sigma), i, width_2;

	width_2 = width / 2;

	den = 2.0*sigma*sigma;
	mask = (float*)calloc(width, sizeof(float));
	for (i = 0; i<width; i++)
	{
		x = (double)(i - width_2);
		mask[i] = (float)exp(-(x*x) / den);
		norm += mask[i];
	}
	for (i = 0; i<width; i++)
		mask[i] /= (float)norm;
  if ( size )
	  *size = width;
	return mask;
}

float *gaussianGradMask1D_float(double sigma, long *size)
{
	float *mask = NULL;
	double norm = 0.0, x, den;
	long width = gaussSigma2Size(sigma), i, width_2;

	width_2 = width / 2;

	den = 2.0*sigma*sigma;
	mask = (float*)calloc(width, sizeof(float));
	for (i = 0; i<width; i++)
	{
		x = (double)(i - width_2);
		mask[i] = (float)(-x*exp(-(x*x) / den));
		// norm += fabs(mask[i]);	
		norm += mask[i] * (width / 2 - i);
	}
	for (i = 0; i<width; i++)
		mask[i] /= (float)norm;
  if ( size )
    *size = width;
	return mask;
}


double *gaussianMask2D(double sigmax, double sigmay, long *sizex, long *sizey)
{
	double *mask = NULL, norm = 0.0, denx, deny, arg;
	long width = gaussSigma2Size(sigmax), 
         height = gaussSigma2Size(sigmay), i, j, width2, height2;

	width2 = width/2;
    height2 = height/2;
    denx = 2.0*sigmax*sigmax;
    deny = 2.0*sigmay*sigmay;
    mask = (double*)calloc(width*height, sizeof(double));
    
    for (j=-height2; j<=height2; j++)
		for (i=-width2; i<=width2; i++)
        {
            arg = (double)(i*i)/denx + (double)(j*j)/deny;
			mask[width*(j+height2)+i+width2] = exp(-arg);
            norm += mask[width*(j+height2)+i+width2];
        }
	for (i=0; i<width*height; i++)
		mask[i] /= norm;
    if ( sizex ) *sizex = width;
    if ( sizey ) *sizey = height;
	return mask;
}

double *gaussianMask3D(double sigmax, double sigmay, double sigmaz, long *sizex, long *sizey, long *sizez)
{
	double *mask = NULL, norm = 0.0, denx, deny, denz, arg;
	long width = gaussSigma2Size(sigmax), 
         height = gaussSigma2Size(sigmay), 
		 depth = gaussSigma2Size(sigmaz),
		 i, j, k, width2, height2, depth2, add;
	
	width2 = width/2;
	height2 = height/2;
	depth2 = depth/2;
	
	denx = 2.0*sigmax*sigmax;
	deny = 2.0*sigmay*sigmay;
	denz = 2.0*sigmaz*sigmaz;
	mask = (double*)calloc(width*height*depth, sizeof(double));
	
	for (k=-depth2; k<=depth2; k++)
		for (j=-height2; j<=height2; j++)
			for (i=-width2; i<=width2; i++)
			{				
				add = ADD3d(i+width2, j+height2, k+depth2, width, height, depth);
				arg = (double)(i*i)/denx + (double)(j*j)/deny + (double)(k*k)/deny;
				mask[add] = exp(-arg);
				norm += mask[add]; 
			}
	
	for (i=0; i<width*height*depth; i++)
		mask[i] /= norm;
	if ( sizex ) *sizex = width;
	if ( sizey ) *sizey = height;
	if ( sizez ) *sizez = depth;
	return mask;
}


void *gaussianMaskFree(void *mask)
{
    FREE(mask)
    return NULL;
}

// ========================================================================

void conv1d(double *data, long size, double *mask, long sizem, double *res)
{
	long i, j, size2;
	double r;

	size2 = sizem/2;

	for (i=0; i<size; i++)
		res[i] = 0.0f;

	for (i=size2; i<size-size2; i++)
	{
		r = 0.0;
		for (j=-size2; j<=size2; j++)
			r += data[i+j] * mask[sizem-(size2+j)-1];	
		res[i] = (double)r;
	}
}

void conv1d_float(float *data, long size, float *mask, long sizem, float *res)
{
	long i, j, size2;
	double r;

	size2 = sizem / 2;

	for (i = 0; i<size; i++)
		res[i] = 0.0f;

	for (i = size2; i<size - size2; i++)
	{
		r = 0.0;
		for (j = -size2; j <= size2; j++)
			r += data[i + j] * mask[sizem - (size2 + j) - 1];
		res[i] = (float)r;
	}
}




void linspace_fill(double x1, double x2, long N, double *v)
{
  long n;
  
  for (n = 0; n < N; n++)
    v[n] = (x2 - x1) / (double)(N - 1)*(double)n + x1;
}


void mesh_fill(double x1, double x2, long Nx,
  double y1, double y2, long Ny,
  double *vx, double *vy)
{
  long x, y;
  double temp;

  for (y=0; y<Ny; y++)
  {
    temp = (y2 - y1) / (double)(Ny - 1)*(double)y + y1;
  
    for (x = 0; x < Nx; x++)
    {
      vx[y*Nx+x] = (x2 - x1) / (double)(Nx - 1)*(double)x + x1;
      vy[y*Nx + x] = temp;
    }
  }
}


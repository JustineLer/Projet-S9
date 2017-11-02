
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <config.h>
#include <eigen.h>
#include <signal.h>
#include <fault_estimation_struct.h>
#include <fault_estimation.h>
#include <util.h>
#include <img_utils.h>
#include <gradient_GPU.cuh>

void data_float_to_out(float * ptr, long width, long height, long depth, char * name);


static long gaussian_size(double sigma)
{
	return 2*(long)ceil(3.0*sigma)+1;
}


float * gaussian_smooth(double sigma, long * size)
{
	float * mask = NULL;
	double v, s = 0.0, x;
	long t, i;

  t = gaussian_size(sigma);
	mask = (float *)calloc(t, sizeof(float));
	for (i=0; i<t; i++)
		{
			x = (double)(i-t/2)/sigma;
			v = exp(-x*x/2.0);
			mask[i] = (float)v;
			s += v;
		}
	for (i=0; i<t; i++) mask[i] /= (float)s;
  if ( size != NULL ) *size = t;
	return mask;
}


float * gaussian_derivative(double sigma, long * size)
{
	float * mask = NULL;
	double v, s = 0.0, x;
	long t, i;

	t = gaussian_size(sigma);
	mask = (float *)calloc(t, sizeof(float));
	for (i=0; i<t; i++)
		{
			x = (double)(i-t/2)/sigma;
			v = -x*exp(-x*x/2.0);
			mask[i] = (float)v;
			s += v*(double)(t/2-i);
		}
	for (i=0; i<t; i++) mask[i] /= (float)s;
  if ( size != NULL ) *size = t;
	return mask;
}


void conv1d(float * v_in, long size, float * mask, long t, float * v_out)
{
	long i, j;
	double s;

	for (i=0; i<size; i++) v_out[i] = 0.0;
	for (i=t/2; i<size-t/2; i++)
		{
			s = 0.0;
			for (j=-t/2; j<=t/2; j++) s += v_in[i+j]*mask[t/2-j];
			v_out[i] = (float)s;
		}
}


static void conv3d_x(float * in, long width, long height, long depth,
							       float * mask, long mask_size, float * v_in, float * v_out,
							       float * out)
{
  long y, z, width2;
  float * v_out2 = NULL;

  width2 = width-mask_size+1;
  v_out2 = v_out+(mask_size-1)/2;
  for (z=0; z<depth; z++)
    for (y=0; y<height; y++)
			{
				data3d_extractx_float(in, width, height, depth, y, z, v_in);
				conv1d(v_in, width, mask, mask_size, v_out);
				data3d_insertx_float(out, width2, height, depth, y, z, v_out2);
			}
}


static void conv3d_y(float * in, long width, long height, long depth,
										 float * mask, long mask_size,
										 float * v_in, float * v_out,
										 float * out)
{
  long x, z, height2;
  float * v_out2 = NULL;

  height2 = height-mask_size+1;
  v_out2 = v_out+(mask_size-1)/2;
  for (z=0; z<depth; z++)
    for (x=0; x<width; x++)
			{
				data3d_extracty_float(in, width, height, depth, x, z, v_in);
				conv1d(v_in, height, mask, mask_size, v_out);
				data3d_inserty_float(out, width, height2, depth, x, z, v_out2);
			}
}


static void conv3d_z(float * in, long width, long height, long depth,
								     float * mask, long mask_size,
										 float * v_in, float * v_out,
										 float * out)
{
  long x, y, depth2;
  float * v_out2 = NULL;

  depth2 = depth-mask_size+1;
  v_out2 = v_out+(mask_size-1)/2;
  for (y=0; y<height; y++)
    for (x=0; x<width; x++)
			{
				data3d_extractz_float(in, width, height, depth, x, y, v_in);
				conv1d(v_in, depth, mask, mask_size, v_out);
				data3d_insertz_float(out, width, height, depth2, x, y, v_out2);
			}
}


static void gradiend3d(float * data, long width, long height, long depth,
											 float * h_lp, float * h_hp, long h_size,
                       float * gx, float * gy, float * gz)
{
	long width2, height2, depth2;
	float * v_in = NULL, * v_out = NULL, * tmp1 = NULL, * tmp2 = NULL;

	v_in = (float*)calloc(width+height+depth, sizeof(float));
	v_out = (float*)calloc(width+height+depth, sizeof(float));
	tmp1 = (float *)calloc(width*height*depth, sizeof(float));
  tmp2 = (float *)calloc(width*height*depth, sizeof(float));

  width2 = width-h_size+1;
  height2 = height-h_size+1;
  depth2 = depth-h_size+1;
  conv3d_z(data, width, height, depth, h_lp, h_size, v_in, v_out, tmp1);
  conv3d_y(tmp1, width, height, depth2, h_lp, h_size, v_in, v_out, tmp2);
  conv3d_x(tmp2, width, height2, depth2, h_hp, h_size, v_in, v_out, gx);
  conv3d_x(tmp1, width, height, depth2, h_lp, h_size, v_in, v_out, tmp2);
  conv3d_y(tmp2, width2, height, depth2, h_hp, h_size, v_in, v_out, gy);
  conv3d_x(data, width, height, depth, h_lp, h_size, v_in, v_out, tmp1);
  conv3d_y(tmp1, width2, height, depth, h_lp, h_size, v_in, v_out, tmp2);
  conv3d_z(tmp2, width2, height2, depth, h_hp, h_size, v_in, v_out, gz);    
	FREE(v_in)
  FREE(v_out)
	FREE(tmp1)
	FREE(tmp2)
}


static void disk_points_compute(double *array_rho, long array_rho_size, double *array_theta, long array_theta_size, double *array_phi, long array_phi_size, double *array_disk, long array_disk_size, 
                                long *Npoint_per_disk, long *Nangle, long ***X, long ***Y, long ***Z, double **ux, double **uy, double **uz, int ***exclusion)
{
  long n, rho_max = 0, irho, N, idx = 0, itheta, iphi, iy, iitheta, iiphi;
  double *x = NULL, *z = NULL, rho, theta, phi, ux0, uy0, uz0, xx, yy, zz, *y = NULL, ux2, uy2, uz2, dot, th;

  for (n = 0; n<array_rho_size; n++)
    rho_max = MAX(rho_max, (long)ceil(array_rho[n]));

  if (array_disk == NULL) array_disk_size = 1;

  *Npoint_per_disk = ( 1 + 2 * (array_rho_size-1) * ((array_rho_size-1) + 1) ) * array_disk_size;
  *Nangle = array_theta_size * array_phi_size;

  CALLOC_TAB((*X), (*Nangle), (*Npoint_per_disk), long)
  CALLOC_TAB((*Y), (*Nangle), (*Npoint_per_disk), long)
  CALLOC_TAB((*Z), (*Nangle), (*Npoint_per_disk), long)

  *ux = (double*)calloc(*Nangle, sizeof(double));
  *uy = (double*)calloc(*Nangle, sizeof(double));
  *uz = (double*)calloc(*Nangle, sizeof(double));

  x = (double*)calloc(*Npoint_per_disk, sizeof(double));
  z = (double*)calloc(*Npoint_per_disk, sizeof(double));
  y = (double*)calloc(*Npoint_per_disk, sizeof(double));
  if ( exclusion )
    CALLOC_VTAB((*exclusion), (*Nangle), (*Nangle), int)
  th = cos(40.0*PI/180.0);

  for (iy=0; iy<array_disk_size; iy++)
    for (irho = 0; irho < array_rho_size; irho++)
    {
      rho = array_rho[irho];
      if (rho == 0) N = 1; else N = 4 * irho; // nombre de points en fonction du rayon
      if (N == 1)
      {
        x[idx] = 0.0;
        if (array_disk) y[idx] = array_disk[iy]; else y[idx] = 0.0;
        z[idx++] = 0.0;
      }
      else
      {
        for (n = 0; n<N; n++)
        {
          theta = 2.0*PI*n / (double)N;
          x[idx] = rho * cos(theta);
          if (array_disk) y[idx] = array_disk[iy]; else y[idx] = 0.0;
          z[idx++] = rho * sin(theta);
        }
      }
    }

  for (itheta = 0; itheta<array_theta_size; itheta++)
    for (iphi = 0; iphi<array_phi_size; iphi++)
    {
      theta = array_theta[itheta];
      phi = array_phi[iphi];
      ux0 = cos(theta) * sin(phi);
      uy0 = -cos(phi);
      uz0 = sin(theta) * sin(phi);
      if (uy0 > 0.0)
      {
        ux0 = -ux0;
        uy0 = -uy0;
        uz0 = -uz0;
      }

      for (n = 0; n<*Npoint_per_disk; n++)
      {
        xx = x[n];
        yy = y[n];
        zz = z[n];
        (*X)[itheta*array_phi_size + iphi][n] = (long)floor(cos(theta) * cos(phi) * xx - cos(theta) * sin(phi) * yy - sin(theta) * zz + .5);
        (*Y)[itheta*array_phi_size + iphi][n] = (long)floor(sin(phi) * xx + cos(phi) * yy + .5);
        (*Z)[itheta*array_phi_size + iphi][n] = (long)floor(sin(theta) * cos(phi) * xx - sin(theta) * sin(phi) * yy + cos(theta) * zz + .5);
      }
      (*ux)[itheta * array_phi_size + iphi] = ux0;
      (*uy)[itheta * array_phi_size + iphi] = uy0;
      (*uz)[itheta * array_phi_size + iphi] = uz0;

      if ( exclusion != NULL )
      {        
        for ( iitheta = 0; iitheta < array_theta_size; iitheta++ )
          for ( iiphi = 0; iiphi < array_phi_size; iiphi++ )
          {
            theta = array_theta[iitheta];
            phi = array_phi[iiphi];
            ux2 = cos(theta) * sin(phi);
            uy2 = -cos(phi);
            uz2 = sin(theta) * sin(phi);
            if ( uy2 > 0.0 )
            {
              ux2 = -ux2;
              uy2 = -uy2;
              uz2 = -uz2;
            }
            dot = ux0*ux2 + uy0*uy2 + uz0*uz2;
            if ( dot > th )
              (*exclusion)[itheta*array_phi_size+iphi][iitheta*array_phi_size + iiphi] = 1;
          }
      }
    }

  FREE(x)
  FREE(y)
  FREE(z)
}

static void best_normal_find(double ux1, double uy1, double uz1, double *ux0, double *uy0, double *uz0, long size, long *nm)
{
  long n;
  double dmax = 0.0, d;

  *nm = 0;
  for ( n = 0; n<size; n++ )
  {
    d = ux1 * ux0[n] + uy1 * uy0[n] + uz1 * uz0[n];
    if ( n == 0 || d > dmax )
    {
      dmax = d;
      *nm = n;
    }
  }
}


static FAULT_ESTIMATION_PARAM *fault_estimation_param3d_init(FAULT_ESTIMATION *fault)
{
  FAULT_ESTIMATION_PARAM *param = NULL;
	FAULT_ESTIMATION_BLOCK_PARAM * block;
  long n, size, border_search = 0;
  
  param = (FAULT_ESTIMATION_PARAM*)calloc(1, sizeof(FAULT_ESTIMATION_PARAM));
	block = &(param->block);
  if ( fault->array_ortho1 )
    {
      for ( n = 0; n < fault->array_ortho1_size; n++ )
        border_search = MAX(border_search, (long)ceil(fabs(fault->array_ortho1[n])));
    }
  block->border1 = (long)ceil(fault->array_rho1[fault->array_rho1_size - 1]) + border_search;
  block->border2 = (long)ceil(fault->array_rho0[fault->array_rho0_size - 1]);
  block->border3 = (long)ceil(3.0*fault->sigma_grad);

  block->mask_grad_lp = gaussian_smooth(fault->sigma_grad, &param->block.mask_grad_size);
  block->mask_grad_hp = gaussian_derivative(fault->sigma_grad, NULL);    

  for ( n = 0; n < 3; n++ )
    block->scale0_size[n] = fault->block_size[n] + 2 * block->border1;
  for ( n = 0; n < 3; n++ )
    block->grad_size[n] = block->scale0_size[2] + 2 * block->border2;
  param->grad_border = (long)ceil(3.0*fault->sigma_grad);
  for ( n = 0; n < 3; n++ )
    block->block_size0[n] = block->grad_size[n] + 2 * block->border3;

  size = fault->block_size[0] * fault->block_size[1] * fault->block_size[2];
  block->crit = (float*)calloc(size, sizeof(float));
  block->ux = (float*)calloc(size, sizeof(float));
  block->uy = (float*)calloc(size, sizeof(float));
  block->uz = (float*)calloc(size, sizeof(float));

  // TODO
  block->crit2 = (float*)calloc(size, sizeof(float));
  block->ux2 = (float*)calloc(size, sizeof(float));
  block->uy2 = (float*)calloc(size, sizeof(float));
  block->uz2 = (float*)calloc(size, sizeof(float));
    
  size = block->scale0_size[0] * block->scale0_size[1] * block->scale0_size[2];
  CALLOC_VTAB(block->scale0_crit, fault->array_phi1_size*fault->array_theta1_size, size, double)      
  
  // disks
  disk_points_compute(fault->array_rho0,  fault->array_rho0_size, fault->array_theta0, fault->array_theta0_size, fault->array_phi0, fault->array_phi0_size, fault->array_ortho0, fault->array_ortho0_size,
                      &param->Npoint_per_disk0, &param->Nangle0, &param->X0, &param->Y0, &param->Z0,
                      &param->ux0, &param->uy0, &param->uz0, NULL);

  disk_points_compute(fault->array_rho1, fault->array_rho1_size, fault->array_theta1, fault->array_theta1_size, fault->array_phi1, fault->array_phi1_size, NULL, 1,
    &param->Npoint_per_disk1, &param->Nangle1, &param->X1, &param->Y1, &param->Z1, &param->ux1, &param->uy1, &param->uz1, &param->array_exclusion);
  
  param->scale1_to_scale0_angle_match = (long*)calloc(param->Nangle1, sizeof(long));
  for (n=0; n<param->Nangle1; n++)
    best_normal_find(param->ux1[n], param->uy1[n], param->uz1[n], param->ux0, param->uy0, param->uz0, param->Nangle0, &param->scale1_to_scale0_angle_match[n]);
 
  // TODO
  return param;
}


void fault_estimation_set_dims(void *_fault, long width, long height, long depth)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  fault->width = width;
  fault->height = height;
  fault->depth = depth;
}

void fault_estimation_get_dims(void *_fault, long *width, long *height, long *depth)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  if (width) *width = fault->width;
  if (height) *height = fault->height;
  if (depth) *depth = fault->depth;
}


void fault_estimation_set_sigma_grad(void *_fault, double sigma)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  fault->sigma_grad = sigma;
}


void fault_estimation_set_block_size(void *_fault, long bx, long by, long bz)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  fault->block_size[0] = bx;
  fault->block_size[1] = by;
  fault->block_size[2] = bz;
}


float *fault_estimation_get_data_out(void *_fault, long type, long *width, long *height, long *depth)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL || fault->param == NULL) return NULL;

  if ( width ) *width = fault->block_size[0];
  if ( height ) *height = fault->block_size[0];
  if ( depth ) *depth = fault->block_size[0];
  if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NX ) return fault->param->block.ux;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NY ) return fault->param->block.uy;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NZ ) return fault->param->block.uz;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_CRIT ) return fault->param->block.crit;
  if ( width ) *width = 0;
  if ( height ) *height = 0;
  if ( depth ) *depth = 0;
  return NULL;
}


static double tensor_crit3d(double sxx, double sxy, double sxz,
  double syy, double syz,
  double szz)
{
  double vp[3];

  eigenCardan(3, sxx, syy, szz, sxy, sxz, syz, vp, NULL, NULL, NULL);
  return (vp[1] + vp[2]) ;
}

static double pca_0(double cn, float *gx, float *gy, float *gz, long width, long height, long depth,
                    long x, long y, long z,
                    long *xx0, long *yy0, long *zz0, long size)
{
  long m, xf, yf, zf, add;
  double sxx = 0.0, sxy = 0.0, sxz = 0.0, syy = 0.0, syz = 0.0, szz = 0.0, c0;

  for (m = 0; m < size; m++)
  {
    xf = xx0[m] + x;
    yf = yy0[m] + y;
    zf = zz0[m] + z;
    if (xf >= 0 && xf < width && yf >= 0 && yf < height && zf >= 0 && zf < depth)
    {
      add = width*height*zf + width*yf + xf;
      sxx += gx[add] * gx[add];
      sxy += gx[add] * gy[add];
      sxz += gx[add] * gz[add];
      syy += gy[add] * gy[add];
      syz += gy[add] * gz[add];
      szz += gz[add] * gz[add];
     }
   }
   sxx /= (cn*size);
   sxy /= (cn*size);
   sxz /= (cn*size);
   syy /= (cn*size);
   syz /= (cn*size);
   szz /= (cn*size);
   c0 = tensor_crit3d(sxx, sxy, sxz, syy, syz, szz);
   return c0;
}


static void fault_block_detection(FAULT_ESTIMATION *fault, long x, long y, long z,
                                  double *c, double *ux, double *uy, double *uz, int *array_exclusion, long *nmax)
{
  long n, nm, nn, *xx0 = NULL, *yy0 = NULL, *zz0 = NULL, offset;
  double ux1, uy1, uz1, cc, cm, c0;
  FAULT_ESTIMATION_PARAM * param = fault->param;

  offset = param->block.border2;
  cm = -1.0;
  *ux = 1.0;
  *uy = 0.0;
  *uz = 0.0;
  if ( nmax ) *nmax = 0;
  for ( n = 0; n < param->Nangle1; n++ )
  {
    if ( array_exclusion == NULL || array_exclusion[n] == 0 )
    {
      ux1 = param->ux1[n];
      uy1 = param->uy1[n];
      uz1 = param->uz1[n];
      nm = param->scale1_to_scale0_angle_match[n];
      xx0 = param->X0[nm];
      yy0 = param->Y0[nm];
      zz0 = param->Z0[nm];
      cc = 1.0;
      for ( nn = 0; nn < param->Npoint_per_disk1; nn++ )
      {
        c0 = fault->param->block.scale0_crit[nm][param->block.scale0_size[0] * param->block.scale0_size[1] * (z + param->Z1[n][nn]) + param->block.scale0_size[0] * (y + param->Y1[n][nn]) + x + param->X1[n][nn]];
        cc *= c0;
      }
      if ( cc > cm )
      {
        cm = cc;
        *ux = ux1;
        *uy = uy1;
        *uz = uz1;
        if ( nmax ) *nmax = n;
      }
    }
  }
  cm = MAX(cm, 0.0);
  *c = pow(cm, 1.0 / param->Npoint_per_disk1);
}


static void fault_estimation_block_scale1_run(void *_fault)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  FAULT_ESTIMATION_PARAM *param = NULL;
  long x, y, z, offset, n, s, add, nmax, nm;
  double ux0, uy0, uz0, critm, uxm, uym, uzm, ux, uy, uz, crit;

  if ( fault == NULL ) return;
  param = fault->param;
  if ( param == NULL ) return;
  
  offset = param->block.border1;
  for ( z = 0; z < fault->block_size[2]; z++ )
    for ( y = 0; y < fault->block_size[1]; y++ )
      for ( x = 0; x < fault->block_size[0]; x++ )
      {
        
        fault_block_detection(fault, x + offset, y + offset, z + offset, &critm, &ux0, &uy0, &uz0, NULL, &nm);
        uxm = ux0;
        uym = uy0;
        uzm = uz0;
        nmax = nm;
        if ( fault->array_ortho1 )
          for ( n = 0; n < fault->array_ortho1_size; n++ )
          {
            s = (long)floor(fault->array_ortho1[n]);
            if ( s == 0.0 ) continue;
            fault_block_detection(fault, (long)(x + offset + s*ux0), (long)(y + offset + s*uy0), (long)(z + offset + s*uz0), &crit, &ux, &uy, &uz, NULL, &nm);
            if ( crit > critm )
            {
              critm = crit;
              uxm = ux;
              uym = uy;
              uzm = uz;
              nmax = nm;
            }
          }
        add = fault->block_size[0] * fault->block_size[1] * z + fault->block_size[0] * y + x;
        fault->param->block.ux[add] = (float)uxm;
        fault->param->block.uy[add] = (float)uym;
        fault->param->block.uz[add] = (float)uzm;
        fault->param->block.crit[add] = (float)critm;
      }
}


static void fault_estimation_scale0_run(void * _fault, float * gx, float * gy, float * gz)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  FAULT_ESTIMATION_PARAM *param = NULL;
  long x, y, z, n, offset;
  double c;

  if ( fault == NULL ) return;
  param = fault->param;
  if ( param == NULL ) return;
  offset = fault->param->block.border2;

  for ( z = 0; z < fault->param->block.scale0_size[2]; z++ )
    for ( y = 0; y < param->block.scale0_size[1]; y++ )
      for ( x = 0; x < param->block.scale0_size[0]; x++ )
        for ( n = 0; n < fault->array_phi0_size*fault->array_theta0_size; n++ )
        {
          c = pca_0(fault->crit_constant+1., gx, gy, gz, param->block.grad_size[0], param->block.grad_size[1], param->block.grad_size[2], x + offset, y + offset, z + offset,
                    param->X0[n], param->Y0[n], param->Z0[n], param->Npoint_per_disk0);
          param->block.scale0_crit[n][param->block.scale0_size[0] * param->block.scale0_size[1] * z + param->block.scale0_size[0] * y + x] = c;
        }
}


void data_crop(float *in, int width_in, int height_in, int depth_in, int b, float *out)
{
	int width_out = width_in - 2 * b, height_out = height_in - 2*b, depth_out = depth_in-2*b, x, y, z, xx, yy, zz;

	for (z = 0; z < depth_out; z++)
	{
		zz = z + b;
		for (y = 0; y < height_out; y++)
		{
			yy = y + b;
			for (x = 0; x < width_out; x++)
			{
				xx = x + b;
				out[width_out*height_out*z + width_out*y + x] = in[width_in*height_in*zz + width_in*yy + xx];
			}
		}
	}
}

void fault_estimation_run(void * fault_estimation, float * data)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
  FAULT_ESTIMATION_PARAM * param = NULL;
  FAULT_ESTIMATION_BLOCK_PARAM * block = NULL;
  long size, x, y, z, xx, yy, zz, b, w, h, d;
  float * tmp = NULL, * gx = NULL, * gy = NULL, * gz = NULL, *gx_t = NULL, *gy_t = NULL, *gz_t = NULL;

  if ( fault_estimation == NULL ) return;
  if ( _fault_estimation->param == NULL )
    _fault_estimation->param = fault_estimation_param3d_init(_fault_estimation);
  if ( data == NULL ) return;
  param = _fault_estimation->param;
  block = &(param->block);
  if ( block == NULL ) return;

	b = block->border1+block->border2+block->border3;
  size = block->block_size0[0]*block->block_size0[1]*block->block_size0[2];
	tmp = (float *)calloc(size, sizeof(float));
	
	w = _fault_estimation->width;
	h = _fault_estimation->height;
	d = _fault_estimation->depth;
  for (z=0; z<block->block_size0[2]; z++)
		{
			if ( z < b ) zz = b-z-1; else if ( z > d-1+b ) zz = 2*d-1+b-z; else zz = z-b;
			for (y = 0; y<block->block_size0[1]; y++)
				{
					if ( y < b ) yy = b-y-1; else if ( y > h-1+b ) yy = 2*h-1+b-y; else yy = y-b;
					for (x=0; x<block->block_size0[0]; x++)
						{
							if ( x < b ) xx = b-x-1; else if ( x > w-1+b ) xx = 2*w-1+b-x; else xx = x-b;
							tmp[block->block_size0[0]*block->block_size0[1]*z+block->block_size0[0]*y+x] = data[w*h*zz+w*yy+xx];
						}			    
				}
		}
  

  gx_t = (float*)calloc(size, sizeof(float));
  gy_t = (float*)calloc(size, sizeof(float));
  gz_t = (float*)calloc(size, sizeof(float));

   gradient_gpu_main(tmp, block->block_size0[0], block->block_size0[1], block->block_size0[2],
	  block->mask_grad_lp, block->mask_grad_hp, block->mask_grad_size, gx_t, gy_t, gz_t);
   size = block->grad_size[0] * block->grad_size[1] * block->grad_size[2];
   gx = (float*)calloc(size, sizeof(float));
   gy = (float*)calloc(size, sizeof(float));
   gz = (float*)calloc(size, sizeof(float));
   b = (block->block_size0[0] - block->grad_size[0]) / 2;
   data_crop(gx_t, block->block_size0[0], block->block_size0[1], block->block_size0[2], b, gx);
   data_crop(gy_t, block->block_size0[0], block->block_size0[1], block->block_size0[2], b, gy);
   data_crop(gz_t, block->block_size0[0], block->block_size0[1], block->block_size0[2], b, gz);
   FREE(gx_t)
   FREE(gy_t)
   FREE(gz_t)

  //gradiend3d(tmp, block->block_size0[0], block->block_size0[1], block->block_size0[2],
             //block->mask_grad_lp, block->mask_grad_hp, block->mask_grad_size, gx, gy, gz);
  
  fault_estimation_scale0_run(_fault_estimation, gx, gy, gz);
	FREE(gx)
	FREE(gy)
	FREE(gz)

  fault_estimation_block_scale1_run(_fault_estimation);
	FREE(tmp)
}


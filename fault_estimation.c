
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


static float compute_max_absolute(float * v, long size)
{
	float v_max = (float)fabs(v[0]);
	long i;

	for (i=1; i<size; i++)
		if ( (float)fabs(v[i]) > v_max ) v_max = (float)fabs(v[i]);
	return v_max;
}


static void disk_points_compute(float * array_rho, long array_rho_size,
																float * array_theta, long array_theta_size,
																float * array_phi, long array_phi_size,
																float * array_disk, long array_disk_size, 
                                long * nb_points, long * nb_angle,
																long *** X, long *** Y, long *** Z,
																float ** ux, float ** uy, float ** uz)
{
  long n, rho_max = 0, irho, nb, idx = 0, itheta, iphi, iy, add;
  float * x = NULL, * y = NULL, * z = NULL, rho, theta, phi, xx, yy, zz, c_theta, s_theta, c_phi, s_phi;

  rho_max = (long)ceil(compute_max_absolute(array_rho, array_rho_size));
  
  *nb_points = (1+2*(array_rho_size-1)*((array_rho_size-1)+1))*array_disk_size;
  *nb_angle = array_theta_size*array_phi_size;

  CALLOC_TAB((*X), (*nb_angle), (*nb_points), long)
  CALLOC_TAB((*Y), (*nb_angle), (*nb_points), long)
  CALLOC_TAB((*Z), (*nb_angle), (*nb_points), long)

  *ux = (float*)calloc(*nb_angle, sizeof(float));
  *uy = (float*)calloc(*nb_angle, sizeof(float));
  *uz = (float*)calloc(*nb_angle, sizeof(float));

  x = (float*)calloc(*nb_points, sizeof(float));
  z = (float*)calloc(*nb_points, sizeof(float));
  y = (float*)calloc(*nb_points, sizeof(float));
  
  for (iy=0; iy<array_disk_size; iy++)
    for (irho=0; irho<array_rho_size; irho++)
			{
				rho = array_rho[irho];
				if (rho == 0) nb = 1; else nb = 4*irho; // nombre de points en fonction du rayon
				if ( nb == 1)
					{
						x[idx] = 0.0f;
						if ( array_disk != NULL ) y[idx] = array_disk[iy]; else y[idx] = 0.0f;
						z[idx++] = 0.0f;
					}
				else
					{
						for (n=0; n<nb; n++)
							{
								theta = (float)(2.0*PI)*(float)n/(float)nb;
								x[idx] = rho*(float)cos((double)theta);
								if ( array_disk != NULL ) y[idx] = array_disk[iy]; else y[idx] = 0.0f;
								z[idx++] = rho*(float)sin(theta);
							}
					}
			}

  for (itheta=0; itheta<array_theta_size; itheta++)
    for (iphi=0; iphi<array_phi_size; iphi++)
			{
				theta = array_theta[itheta];
				c_theta = (float)cos(theta);
				s_theta = (float)sin(theta);
				phi = array_phi[iphi];
				c_phi = (float)cos(phi);
				s_phi = (float)sin(phi);
				add = itheta*array_phi_size+iphi;
				(*ux)[add] = c_theta*s_phi;
				(*uy)[add] = -c_phi;
				(*uz)[add] = s_theta*s_phi;
				for (n=0; n<*nb_points; n++)
					{
						xx = x[n];
						yy = y[n];
						zz = z[n];
						(*X)[add][n] = (long)floor((double)(c_theta*c_phi*xx-c_theta*s_phi*yy-s_theta*zz+0.5f));
						(*Y)[add][n] = (long)floor((double)(s_phi*xx+c_phi*yy+0.5f));
						(*Z)[add][n] = (long)floor((double)(s_theta*c_phi*xx-s_theta*s_phi*yy+c_theta*zz+0.5f));
					}
			}
  FREE(x)
  FREE(y)
  FREE(z)
}


static float st(float * gx, float * gy, float * gz, long sx, long sy, long sz,
								long x, long y, long z, long * dx, long * dy, long * dz, long size)
{
  long n, xf, yf, zf, add;
  double sxx = 0.0, sxy = 0.0, sxz = 0.0, syy = 0.0, syz = 0.0, szz = 0.0, vp[3];

  for (n=0; n<size; n++)
		{
			xf = dx[n]+x;
			yf = dy[n]+y;
			zf = dz[n]+z;
			add = sx*sy*zf+sx*yf+xf;
			sxx += (float)(gx[add]*gx[add]);
			sxy += (float)(gx[add]*gy[add]);
			sxz += (float)(gx[add]*gz[add]);
			syy += (float)(gy[add]*gy[add]);
			syz += (float)(gy[add]*gz[add]);
			szz += (float)(gz[add]*gz[add]);
		}
  eigenCardan(3, sxx, syy, szz, sxy, sxz, syz, vp, NULL, NULL, NULL);
  return (float)((vp[1]+vp[2])/(double)size);
}


static long crit_combine(float ** v_crit, long sx, long sy, long sz,
												 long * scale1_to_scale0,
												 long x, long y, long z,
												 long ** dx, long ** dy, long ** dz, long nb_angle, long nb_points,
												 float * crit_opt)
{
  long n, p, first = 1, n_max = 0;
  float crit, crit_max, * v_crit_n = NULL;
	
  for (n=0; n<nb_angle; n++)
		{
			v_crit_n = v_crit[scale1_to_scale0[n]];
			crit = 1.0;
			for (p=0; p<nb_points; p++)
				crit *= v_crit_n[sx*sy*(z+dz[n][p])+sx*(y+dy[n][p])+x+dx[n][p]];
			if ( first || crit > crit_max )
				{
					first = 0;
					crit_max = crit;
					n_max = n;
				}
		}
	*crit_opt = (float)pow((double)crit_max, 1.0/(double)nb_points);
	return n_max;
}


static void scale1_run(float ** v_crit, long nb_angle0, long sx0, long sy0, long sz0,
											 long * scale1_to_scale0, long offset,
											 float * ux1, float * uy1, float * uz1,
											 long ** dx, long ** dy, long ** dz, long nb_angle1, long nb_points1,
											 float * array_ortho, long array_ortho_size,
											 float * ux, float * uy, float * uz, float * crit,
											 long sx1, long sy1, long sz1)
{
  long x, y, z, n, add, n_c, n_e, xx, yy, zz;
  float crit_c, crit_e;

  for (z=0; z<sz1; z++)
    for (y=0; y<sy1; y++)
      for (x=0; x<sx1; x++)
				{
					n_c = crit_combine(v_crit, sx0, sy0, sz0, scale1_to_scale0, x+offset, y+offset, z+offset,
														 dx, dy, dz, nb_angle1, nb_points1, &crit_c);
					if ( array_ortho )
						for (n=0; n<array_ortho_size; n++)
							{
								xx = (long)floor((double)(array_ortho[n]*ux1[n_c])+0.5);
								yy = (long)floor((double)(array_ortho[n]*uy1[n_c])+0.5);
								zz = (long)floor((double)(array_ortho[n]*uz1[n_c])+0.5);
								n_e = crit_combine(v_crit, sx0, sy0, sz0, scale1_to_scale0, x+xx+offset, y+yy+offset, z+zz+offset,
																	 dx, dy, dz, nb_angle1, nb_points1, &crit_e);
								if ( crit_e > crit_c )
									{
										crit_c = crit_e;
										n_c = n_e;
									}
							}
					add = sx1*sy1*z+sx1*y+x;
					ux[add] = ux1[n_c];
					uy[add] = uy1[n_c];
					uz[add] = uz1[n_c];
					crit[add] = crit_c;
				}
}


static void scale0_run(float * gx, float * gy, float * gz, long offset,
											 long ** dx, long ** dy, long ** dz, long nb_angle, long nb_points,
											 float ** v_crit, long sx, long sy, long sz)
{
	long x, y, z, n;

  for (z=0; z<sz; z++)
    for (y=0; y<sy; y++)
      for (x=0; x<sx; x++)
        for (n=0; n<nb_angle; n++)
					v_crit[n][sx*sy*z+sx*y+x]
						= st(gx, gy, gz, sx+offset*2, sy+offset*2, sz+offset*2, x+offset, y+offset, z+offset, dx[n], dy[n], dz[n], nb_points);
}


static long close_normal(float ux, float uy, float uz, float * vx, float * vy, float * vz, long size)
{
  double d, d_max;
	long n, n_max;

  for (n=0; n<size; n++)
		{
			d = ux*vx[n]+uy*vy[n]+uz*vz[n];
			if ( n == 0 || d > d_max )
				{
					d_max = d;
					n_max = n;
				}
		}
	return n_max;
}


void fault_estimation_get_borders(void * fault_estimation, long * v_border)
{
	FAULT_ESTIMATION * f = (FAULT_ESTIMATION *)fault_estimation;

  if ( fault_estimation == NULL ) return;
	v_border[0] = gaussian_size(f->sigma_grad)/2;
	v_border[1] = (long)ceil((double)compute_max_absolute(f->array_ortho0, f->array_ortho0_size));
	v_border[1] += (long)ceil((double)compute_max_absolute(f->array_rho0, f->array_rho0_size));
	v_border[2] = (long)ceil((double)compute_max_absolute(f->array_ortho1, f->array_ortho1_size));
	v_border[2] += (long)ceil((double)compute_max_absolute(f->array_rho1, f->array_rho1_size));
 }


static void fault_estimation_init(FAULT_ESTIMATION * fault_estimation)
{
  FAULT_ESTIMATION * f = (FAULT_ESTIMATION *)fault_estimation;
	long n, size;
	float * ux0 = NULL, * uy0 = NULL, * uz0 = NULL;
  
	size = f->block_size[0]*f->block_size[1]*f->block_size[2];
  f->crit = (float *)calloc(size, sizeof(float));
  f->nx = (float *)calloc(size, sizeof(float));
  f->ny = (float *)calloc(size, sizeof(float));
  f->nz = (float *)calloc(size, sizeof(float));
  
  disk_points_compute(f->array_rho0,  f->array_rho0_size,
                      f->array_theta0, f->array_theta0_size,
                      f->array_phi0, f->array_phi0_size,
                      f->array_ortho0, f->array_ortho0_size,
                      &f->nb_points0, &f->nb_angle0,
                      &f->X0, &f->Y0, &f->Z0,
                      &ux0, &uy0, &uz0);

  disk_points_compute(f->array_rho1, f->array_rho1_size,
											f->array_theta1, f->array_theta1_size,
											f->array_phi1, f->array_phi1_size,
											NULL, 1,
											&f->nb_points1, &f->nb_angle1,
											&f->X1, &f->Y1, &f->Z1,
											&f->ux1, &f->uy1, &f->uz1);
  
  f->scale1_to_scale0 = (long *)calloc(f->nb_angle1, sizeof(long));
  for (n=0; n<f->nb_angle1; n++)
    f->scale1_to_scale0[n] = close_normal(f->ux1[n], f->uy1[n], f->uz1[n], ux0, uy0, uz0, f->nb_angle0);

	FREE(ux0)
  FREE(uy0)
  FREE(uz0)
}


void fault_estimation_run(void * fault_estimation, float * data_in)
{
  FAULT_ESTIMATION * f = (FAULT_ESTIMATION *)fault_estimation;
  long v_border[3], w, h, d, size, grad_size, sx0, sy0, sz0, size0, sx1, sy1, sz1;
  float * tmp = NULL, * gx = NULL, * gy = NULL, * gz = NULL, * grad_lp = NULL, * grad_hp = NULL, ** scale0_crit = NULL;
  
	if ( data_in == NULL ) return;
  if ( fault_estimation == NULL ) return;
	if ( f->init == 0 )
		{
			f->init = 1;
			fault_estimation_init(f);
		}
   
	fault_estimation_get_borders(fault_estimation, v_border);
  
	w = f->block_size[0]+2*(v_border[0]+v_border[1]+v_border[2]);
	h = f->block_size[1]+2*(v_border[0]+v_border[1]+v_border[2]);
	d = f->block_size[2]+2*(v_border[0]+v_border[1]+v_border[2]);
	size = (w-2*v_border[0])*(h-2*v_border[0])*(d-2*v_border[0]);

	gx = (float *)calloc(size, sizeof(float));
	gy = (float *)calloc(size, sizeof(float));
	gz = (float *)calloc(size, sizeof(float));

	grad_lp = gaussian_smooth(f->sigma_grad, &grad_size);
	grad_hp = gaussian_derivative(f->sigma_grad, NULL); 
	if ( f->cpugpu == 0 )  
		gradiend3d(data_in, w, h, d, grad_lp, grad_hp, grad_size, gx, gy, gz);
	else
	{
		gradient_gpu_main(data_in, w, h, d, grad_lp, grad_hp, grad_size, gx, gy, gz);
		data_float_to_out(gx, w - 2 * v_border[0], h - 2 * v_border[0], d - 2 * v_border[0], "gx");
		data_float_to_out(gy, w - 2 * v_border[0], h - 2 * v_border[0], d - 2 * v_border[0], "gy");
		data_float_to_out(gz, w - 2 * v_border[0], h - 2 * v_border[0], d - 2 * v_border[0], "gz");
	}

	FREE(grad_lp)
	FREE(grad_hp)

	sx0 = f->block_size[0]+2*v_border[2];
	sy0 = f->block_size[1]+2*v_border[2];
	sz0 = f->block_size[2]+2*v_border[2];
	size0 = sx0*sy0*sz0;
	CALLOC_VTAB(scale0_crit, f->array_phi0_size*f->array_theta0_size, size0, float) 
	if ( f->cpugpu == 0 )
		scale0_run(gx, gy, gz, v_border[1], f->X0, f->Y0, f->Z0, f->array_phi0_size*f->array_theta0_size, f->nb_points0,
							 scale0_crit, sx0, sy0, sz0);
	else 
		scale0_run(gx, gy, gz, v_border[1], f->X0, f->Y0, f->Z0, f->array_phi0_size*f->array_theta0_size, f->nb_points0,
			scale0_crit, sx0, sy0, sz0);// A faire en version GPU
	FREE(gx)
	FREE(gy)
	FREE(gz)

	sx1 = f->block_size[0];
	sy1 = f->block_size[1];
	sz1 = f->block_size[2];
	if ( f->cpugpu == 0 )
		scale1_run(scale0_crit, f->array_phi0_size*f->array_theta0_size, sx0, sy0, sz0, f->scale1_to_scale0, v_border[2],
							 f->ux1, f->uy1, f->uz1,
							 f->X1, f->Y1, f->Z1, f->array_phi1_size*f->array_theta1_size, f->nb_points1,
							 f->array_ortho1, f->array_ortho1_size,	
							 f->nx, f->ny, f->nz, f->crit,
							 sx1, sy1, sz1);
	else 
		scale1_run(scale0_crit, f->array_phi0_size*f->array_theta0_size, sx0, sy0, sz0, f->scale1_to_scale0, v_border[2],
			f->ux1, f->uy1, f->uz1,
			f->X1, f->Y1, f->Z1, f->array_phi1_size*f->array_theta1_size, f->nb_points1,
			f->array_ortho1, f->array_ortho1_size,
			f->nx, f->ny, f->nz, f->crit,
			sx1, sy1, sz1);// A faire en version GPU


	FREE_VTAB(scale0_crit)
}


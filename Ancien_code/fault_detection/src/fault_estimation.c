
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
// #include <time.h>

#include <util.h>
#include <eigen.h>
#include <chronos.h>
#include <signal.h>
#include <convolution_separable_3d.h>
#include <fault_estimation.h>

static double pca_0(double cn, float *gx, float *gy, float *gz, long width, long height, long depth,
                    long x, long y, long z,
                    long *xx0, long *yy0, long *zz0, long size);


typedef struct _FAULT_ESTIMATION_BLOCK_PARAM
{
  long block_size0[3],    // dimensions block initial
    Nblocks[3],         // nombre de blocks
    grad_size[3], scale0_size[3], border1, border2, border3, mask_grad_size;
  float *gx, *gy, *gz, *ux, *uy, *uz, *crit, *tmp1, *tmp2, *mask_grad_lp, *mask_grad_hp, *ux2, *uy2, *uz2, *crit2;
  double **scale0_crit;

}FAULT_ESTIMATION_BLOCK_PARAM;

typedef struct _FAULT_ESTIMATION_PARAM
{
  long **X0, **Y0, **Z0, Nangle0, Npoint_per_disk0, Npoint_per_disk1, cpt_temp, cpt_temp2, Nangle1, **X1, **Y1, **Z1,
       plane_border[3],   // plus grand pour la recherche des plans
       grad_border,       // plus grand pour le calcul du gradient
       border[3],         // bordure totale (plan + gradient)
       *scale1_to_scale0_angle_match;

  double *ux0, *uy0, *uz0, *ux1, *uy1, *uz1;
  int **array_exclusion;
  // float *gx, *gy, *gz;
  FAULT_ESTIMATION_BLOCK_PARAM *block;
  void *file_nx, *file_ny, *file_nz, *file_crit, *file_nx2, *file_ny2, *file_nz2, *file_crit2;
}FAULT_ESTIMATION_PARAM;


typedef struct _FAULT_ESTIMATION
{
  long width, height, depth, mix_type, nbthreads, 
       array_theta0_size, array_phi0_size, array_theta1_size, array_phi1_size, array_rho0_size, array_rho1_size, size_in, array_disk0_size,
       block_size[3], run_type, array_search1_size;
  double *array_theta0, *array_phi0, *array_theta1, *array_phi1, *array_rho0, *array_rho1, *Xin, *Yin, *Zin, crit_constant, *array_disk0, sigma_grad, *array_search_1;
  float *nx, *ny, *nz, *crit;
  float *gx, *gy, *gz;
  void *chronos;
  char *data_in_filename, *nx_filename, *ny_filename, *nz_filename, *crit_filename;
  FAULT_ESTIMATION_PARAM *param;
}FAULT_ESTIMATION;


static double fault_estimation_gradient_constant_compute(float *gx, float *gy, float *gz, long size)
{
  long n;
  double s = 0.0;

  if (gx == NULL || gy == NULL) return 0.0;
  if (gz == NULL)
  {
    for (n = 0; n<size; n++)
      s += sqrt(SQR(gx[n]) + SQR(gy[n]));
    return SQR(s / (double)(size));
  }
  else
  {
    for (n = 0; n<size; n++)
      s += sqrt(SQR(gx[n]) + SQR(gy[n]) + SQR(gz[n]));
    return SQR(s / (double)(size));
  }
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
            // dot = fabs(ux0*ux2 + uy0*uy2 + uz0*uz2);
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
  long n, size, border_search = 0;
  
  param = (FAULT_ESTIMATION_PARAM*)calloc(1, sizeof(FAULT_ESTIMATION_PARAM));

  // pour les blocks

  if ( fault->run_type == FAULT_ESTIMATION_RUN_TYPE_BLOCK )
  {
    param->block = (FAULT_ESTIMATION_BLOCK_PARAM*)calloc(1, sizeof(FAULT_ESTIMATION_BLOCK_PARAM));
    if ( fault->array_search_1 )
    {
      for ( n = 0; n < fault->array_search1_size; n++ )
        border_search = MAX(border_search, (long)ceil(fabs(fault->array_search_1[n])));
    }
    param->block->border1 = (long)ceil(fault->array_rho1[fault->array_rho1_size - 1]) + border_search;
    param->block->border2 = (long)ceil(fault->array_rho0[fault->array_rho0_size - 1]);
    param->block->border3 = (long)ceil(3.0*fault->sigma_grad);

    param->block->mask_grad_lp = gaussianMask1D_float(fault->sigma_grad, &param->block->mask_grad_size);
    param->block->mask_grad_hp = gaussianGradMask1D_float(fault->sigma_grad, NULL);    

    for ( n = 0; n < 3; n++ )
      param->block->scale0_size[n] = fault->block_size[n] + 2 * param->block->border1;
    for ( n = 0; n < 3; n++ )
      param->block->grad_size[n] = param->block->scale0_size[2] + 2 * param->block->border2;
    param->grad_border = (long)ceil(3.0*fault->sigma_grad);
    for ( n = 0; n < 3; n++ )
      param->block->block_size0[n] = param->block->grad_size[n] + 2 * param->block->border3;

    size = param->block->block_size0[0] * param->block->block_size0[1] * param->block->block_size0[2];
    param->block->tmp1 = (float*)calloc(size, sizeof(float));
    param->block->tmp2 = (float*)calloc(size, sizeof(float));

    // for (n=0; n<3; n++ ) param->plane_border[n] = fault->array_rho1[fault->array_rho1_size - 1] + fault->array_rho0[fault->array_rho0_size - 1];
    
    param->block->Nblocks[0] = (fault->width - 1) / fault->block_size[0] + 1;
    param->block->Nblocks[1] = (fault->height - 1) / fault->block_size[1] + 1;
    param->block->Nblocks[2] = (fault->depth - 1) / fault->block_size[2] + 1;

    size = param->block->grad_size[0] * param->block->grad_size[1] * param->block->grad_size[2];
    param->block->gx = (float*)calloc(size, sizeof(float));
    param->block->gy = (float*)calloc(size, sizeof(float));
    param->block->gz = (float*)calloc(size, sizeof(float));
    
    size = fault->block_size[0] * fault->block_size[1] * fault->block_size[2];
    param->block->crit = (float*)calloc(size, sizeof(float));
    param->block->ux = (float*)calloc(size, sizeof(float));
    param->block->uy = (float*)calloc(size, sizeof(float));
    param->block->uz = (float*)calloc(size, sizeof(float));

    // TODO
    param->block->crit2 = (float*)calloc(size, sizeof(float));
    param->block->ux2 = (float*)calloc(size, sizeof(float));
    param->block->uy2 = (float*)calloc(size, sizeof(float));
    param->block->uz2 = (float*)calloc(size, sizeof(float));
    

    size = param->block->scale0_size[0] * param->block->scale0_size[1] * param->block->scale0_size[2];
    CALLOC_VTAB(param->block->scale0_crit, fault->array_phi1_size*fault->array_theta1_size, size, double)      
  }
  
  // disks
  disk_points_compute(fault->array_rho0,  fault->array_rho0_size, fault->array_theta0, fault->array_theta0_size, fault->array_phi0, fault->array_phi0_size, fault->array_disk0, fault->array_disk0_size,
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

void fault_estimation_end(void *_fault)
{
  // fault_estimation_files_close(_fault);
}


static FAULT_ESTIMATION_BLOCK_PARAM *fault_estimate_block_param_realease(FAULT_ESTIMATION_BLOCK_PARAM *block)
{
  if ( block == NULL ) return NULL;
  FREE_VTAB(block->scale0_crit)
  FREE(block->gx)
  FREE(block->gy)
  FREE(block->gz)
  FREE(block->ux)
  FREE(block->uy)
  FREE(block->uz)
  FREE(block->crit)

  FREE(block->ux2)
  FREE(block->uy2)
  FREE(block->uz2)
  FREE(block->crit2)

  FREE(block->tmp1)
  FREE(block->tmp2)
  block->mask_grad_lp = gaussianMaskFree(block->mask_grad_lp);
  block->mask_grad_hp = gaussianMaskFree(block->mask_grad_hp);
  FREE_VTAB(block->scale0_crit)
  FREE(block)
  return NULL;
}



static FAULT_ESTIMATION_PARAM *fault_estimation_param_release(FAULT_ESTIMATION_PARAM *param)
{
  if (param == NULL) return NULL;
  
  FREE_TAB(param->X0, param->Nangle0)
  FREE_TAB(param->Y0, param->Nangle0)
  FREE_TAB(param->Z0, param->Nangle0)
  FREE_TAB(param->X1, param->Nangle0)
  FREE_TAB(param->Y1, param->Nangle0)
  FREE_TAB(param->Z1, param->Nangle0)
  FREE(param->scale1_to_scale0_angle_match)
  FREE(param->ux0)
  FREE(param->uy0)
  FREE(param->uz0)
  FREE(param->ux1)
  FREE(param->uy1)
  FREE(param->uz1)
  FREE_VTAB(param->array_exclusion)
  param->block = fault_estimate_block_param_realease(param->block);
  FREE(param)
  return NULL;
}


void *fault_estimation_init()
{
  FAULT_ESTIMATION *data = NULL;

  data = (FAULT_ESTIMATION*)calloc(1, sizeof(FAULT_ESTIMATION));
  fault_estimation_set_nbthreads(data, 1);
  fault_estimation_set_run_type(data, FAULT_ESTIMATION_RUN_TYPE_BLOCK);
  fault_estimation_set_block_size(data, 64, 64, 64);
  fault_estimation_set_sigma_grad(data, 1.0);
  // fault_estimation2_set_crit_constant(data, 1.0);
  // fault_estimation2_set_mix_type(data, FAULT_ESTIMATION_MIX_TYPE_NONE);
  return (void*)data;
}

void *fault_estimation_release(void *_fault)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return NULL;
  fault->param = fault_estimation_param_release(fault->param);
  FREE(fault->array_theta0)
  FREE(fault->array_phi0)
  FREE(fault->array_theta1)
  FREE(fault->array_phi1)
  FREE(fault->array_search_1)
  FREE(fault->data_in_filename)
  FREE(fault->nx_filename)
  FREE(fault->ny_filename)
  FREE(fault->nz_filename)
  FREE(fault->crit_filename)
  fault->chronos = CHRONOS_RELEASE(fault->chronos);
  FREE(fault)
  return NULL;
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

void fault_estimation_set_array_rho0(void *_fault, double *rho, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if (fault == NULL) return;
  FREE(fault->array_rho0)
    fault->array_rho0 = (double*)calloc(size, sizeof(double));
  fault->array_rho0_size = size;
  for (n = 0; n<size; n++)
    fault->array_rho0[n] = rho[n];
}

void fault_estimation_set_array_rho1(void *_fault, double *rho, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if (fault == NULL) return;
  FREE(fault->array_rho1)
    fault->array_rho1 = (double*)calloc(size, sizeof(double));
  fault->array_rho1_size = size;
  for (n = 0; n<size; n++)
    fault->array_rho1[n] = rho[n];
}

void fault_estimation_set_array_theta0(void *_fault, double *theta, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if (fault == NULL) return;
  FREE(fault->array_theta0)
  fault->array_theta0 = (double*)calloc(size, sizeof(double));
  fault->array_theta0_size = size;
  for (n = 0; n<size; n++)
    fault->array_theta0[n] = theta[n];
}

void fault_estimation_set_array_phi0(void *_fault, double *phi, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if (fault == NULL) return;
  FREE(fault->array_phi0)
  fault->array_phi0 = (double*)calloc(size, sizeof(double));
  fault->array_phi0_size = size;
  for (n = 0; n<size; n++)
    fault->array_phi0[n] = phi[n];
}

void fault_estimation_set_array_theta1(void *_fault, double *theta, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if (fault == NULL) return;
  FREE(fault->array_theta1)
    fault->array_theta1 = (double*)calloc(size, sizeof(double));
  fault->array_theta1_size = size;
  for (n = 0; n<size; n++)
    fault->array_theta1[n] = theta[n];
}

void fault_estimation_set_array_phi1(void *_fault, double *phi, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if (fault == NULL) return;
  FREE(fault->array_phi1)
    fault->array_phi1 = (double*)calloc(size, sizeof(double));
  fault->array_phi1_size = size;
  for (n = 0; n<size; n++)
    fault->array_phi1[n] = phi[n];
}

void fault_estimation_set_array_search1(void *_fault, double *pos, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if ( fault == NULL ) return;
  FREE(fault->array_search_1)
  fault->array_search_1 = (double*)calloc(size, sizeof(double));
  fault->array_search1_size = size;
  for ( n = 0; n < size; n++ )
    fault->array_search_1[n] = pos[n];
}

void fault_estimation_set_array_disk0(void *_fault, double *array_disk, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long n;

  if (fault == NULL) return;
  FREE(fault->array_disk0)
  fault->array_disk0 = (double*)calloc(size, sizeof(double));
  for (n = 0; n < size; n++)
    fault->array_disk0[n] = array_disk[n];
  fault->array_disk0_size = size;
}


void fault_estimation_set_domain_in(void *_fault, double *X, double *Y, double *Z, long size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  fault->Xin = X;
  fault->Yin = Y;
  fault->Zin = Z;
  fault->size_in = size;
}

void fault_estimation_get_domain_in(void *_fault, double **X, double **Y, double **Z, long *size)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  if (X) *X = fault->Xin;
  if (Y) *Y = fault->Yin;
  if (Z) *Z = fault->Zin;
  if (size) *size = fault->size_in;
}


void fault_estimation_set_gradient(void *_fault, float *gx, float *gy, float *gz)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  fault->gx = gx;
  fault->gy = gy;
  fault->gz = gz;
}

void fault_estimation_set_sigma_grad(void *_fault, double sigma)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  fault->sigma_grad = sigma;
}

//void fault_estimation_get_gradient(void *_fault, float **gx, float **gy, float **gz)
//{
//  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
//
//  if (fault == NULL) return;
//  if (gx) *gx = fault->gx;
//  if (gy) *gy = fault->gy;
//  if (gz) *gz = fault->gz;
//}
//
//

void fault_estimation_set_out(void *_fault, float *nx, float *ny, float *nz, float *crit)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  fault->nx = nx;
  fault->ny = ny;
  fault->nz = nz;
  fault->crit = crit;
}

void fault_estimation_set_out_filename(void *_fault, char *nx_filename, char *ny_filename, char *nz_filename, char *crit_filename)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  ARRAY_CHAR_COPY(nx_filename, fault->nx_filename)
  ARRAY_CHAR_COPY(ny_filename, fault->ny_filename)
  ARRAY_CHAR_COPY(nz_filename, fault->nz_filename)
  ARRAY_CHAR_COPY(crit_filename, fault->crit_filename)
}

void fault_estimation_set_chronos(void *_fault, long val)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  if (val == 0)
    fault->chronos = CHRONOS_RELEASE(fault->chronos);
  else
    fault->chronos = CHRONOS_INIT;
}

void fault_estimation_chronos_reset(void *_fault, long label)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  CHRONOS_RESET(fault->chronos, label);
}

double fault_estimation_chronos_get(void *_fault, long label)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return 0.0;
  return CHRONOS_GET(fault->chronos, label);
}

void fault_estimation_set_crit_constant(void *_fault, double val)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if (fault == NULL) return;
  fault->crit_constant = val;
}

void fault_estimation_set_nbthreads(void *_fault, long nbthreads)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  fault->nbthreads = nbthreads;
}

void fault_estimation_set_data_in_filename(void *_fault, char *filename)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  ARRAY_CHAR_COPY(filename, fault->data_in_filename)
}

void fault_estimation_set_block_size(void *_fault, long bx, long by, long bz)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  fault->block_size[0] = bx;
  fault->block_size[1] = by;
  fault->block_size[2] = bz;
}

void fault_estimation_set_run_type(void *_fault, long type)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return;
  fault->run_type = type;
}

float *fault_estimation_get_data_out(void *_fault, long type, long *width, long *height, long *depth)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL || fault->param == NULL) return NULL;

  if ( width ) *width = fault->block_size[0];
  if ( height ) *height = fault->block_size[0];
  if ( depth ) *depth = fault->block_size[0];
  if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NX ) return fault->param->block->ux;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NY ) return fault->param->block->uy;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NZ ) return fault->param->block->uz;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_CRIT ) return fault->param->block->crit;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NX2 ) return fault->param->block->ux2;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NY2 ) return fault->param->block->uy2;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_NZ2 ) return fault->param->block->uz2;
  else if ( type == FAULT_ESTIMATION_GET_DATA_OUT_CRIT2 ) return fault->param->block->crit2;
  if ( width ) *width = 0;
  if ( height ) *height = 0;
  if ( depth ) *depth = 0;
  return NULL;
}

long fault_estimation_get_conv_border(void *_fault)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;

  if ( fault == NULL ) return 0;
  if ( fault->param == NULL ) return 0;
  if ( fault->param->block == NULL ) return 0;
  return fault->param->block->border1 + fault->param->block->border2 + fault->param->block->border3;
}

void fault_estimation_get_nread_block(void *_fault, int *dim)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  long i;
  
  if ( fault == NULL || fault->param == NULL ) return;
  for (i=0; i<3; i++ )
    dim[i] = fault->param->block->Nblocks[i];
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

static void fault_detection(FAULT_ESTIMATION *fault, long x, long y, long z, double *c, double *ux, double *uy, double *uz)
{
  long n, nm, nn, *xx0 = NULL, *yy0 = NULL, *zz0 = NULL;
  double ux1, uy1, uz1, cc, cm, c0;
  FAULT_ESTIMATION_PARAM *param = fault->param;
  
  cm = 0.0;
  for (n = 0; n < param->Nangle1; n++)
  {
    ux1 = param->ux1[n];
    uy1 = param->uy1[n];
    uz1 = param->uz1[n];
    // best_normal_find(ux1, uy1, uz1, param->ux0, param->uy0, param->uz0, param->Nangle0, &nm);
    nm = param->scale1_to_scale0_angle_match[n];

    xx0 = param->X0[nm];
    yy0 = param->Y0[nm];
    zz0 = param->Z0[nm];

    cc = 1.0;
    for (nn=0; nn<param->Npoint_per_disk1; nn++)
    { 
       c0 = pca_0(fault->crit_constant+1.0, fault->gx, fault->gy, fault->gz, fault->width, fault->height, fault->depth,
                  x+param->X1[n][nn], y+param->Y1[n][nn], z+ param->Z1[n][nn], xx0, yy0, zz0, param->Npoint_per_disk0);
       cc *= c0;                
    }
    if ( cc > cm )
    {
      cm = cc;
      *ux = ux1;
      *uy = uy1;
      *uz = uz1;    
    }
  }
  *c = pow(cm, 1.0/ param->Npoint_per_disk1);
}


static void fault_block_detection(FAULT_ESTIMATION *fault, long x, long y, long z,
                                  double *c, double *ux, double *uy, double *uz, int *array_exclusion, long *nmax)
{
  long n, nm, nn, *xx0 = NULL, *yy0 = NULL, *zz0 = NULL, offset;
  double ux1, uy1, uz1, /*uxm, uym, uzm,*/ cc, cm, c0;
  FAULT_ESTIMATION_PARAM *param = fault->param;

  offset = fault->param->block->border2;
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
      // uxm = param->ux0[nm];
      // uym = param->uy0[nm];
      // uzm = param->uz0[nm];
      xx0 = param->X0[nm];
      yy0 = param->Y0[nm];
      zz0 = param->Z0[nm];
      cc = 1.0;
      for ( nn = 0; nn < param->Npoint_per_disk1; nn++ )
      {
        c0 = fault->param->block->scale0_crit[nm][param->block->scale0_size[0] * param->block->scale0_size[1] * (z + param->Z1[n][nn]) + param->block->scale0_size[0] * (y + param->Y1[n][nn]) + x + param->X1[n][nn]];
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

typedef struct _FAULT_THREAD_DATA
{
  void *_fault;
  long offset, size, no;
  double *Xin, *Yin, *Zin;
}FAULT_THREAD_DATA;



static void fault_estimation_block_scale1_run(void *_fault)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  FAULT_ESTIMATION_PARAM *param = NULL;
  long x, y, z, offset, n, s, add, nmax, nm;
  double ux0, uy0, uz0, critm, uxm, uym, uzm, ux, uy, uz, crit;

  if ( fault == NULL ) return;
  param = fault->param;
  if ( param == NULL ) return;
  
  offset = fault->param->block->border1;
  for ( z = 0; z < fault->block_size[2]; z++ )
    for ( y = 0; y < fault->block_size[1]; y++ )
      for ( x = 0; x < fault->block_size[0]; x++ )
      {
        
        fault_block_detection(fault, x + offset, y + offset, z + offset, &critm, &ux0, &uy0, &uz0, NULL, &nm);
        uxm = ux0;
        uym = uy0;
        uzm = uz0;
        nmax = nm;
        if ( fault->array_search_1 )
          for ( n = 0; n < fault->array_search1_size; n++ )
          {
            s = (long)floor(fault->array_search_1[n]);
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
        fault->param->block->ux[add] = (float)uxm;
        fault->param->block->uy[add] = (float)uym;
        fault->param->block->uz[add] = (float)uzm;
        fault->param->block->crit[add] = (float)critm;

        /*
        fault_block_detection(fault, x + offset, y + offset, z + offset, &critm, &ux0, &uy0, &uz0, param->array_exclusion[nmax], NULL);
        uxm = ux0;
        uym = uy0;
        uzm = uz0;
        if ( fault->array_search_1 )
          for ( n = 0; n < fault->array_search1_size; n++ )
          {
            s = (long)floor(fault->array_search_1[n]);
            if ( s == 0.0 ) continue;
            fault_block_detection(fault, (long)(x + offset + s*ux0), (long)(y + offset + s*uy0), (long)(z + offset + s*uz0), &crit, &ux, &uy, &uz, param->array_exclusion[nmax], NULL);
            if ( crit > critm )
            {
              critm = crit;
              uxm = ux;
              uym = uy;
              uzm = uz;
            }
          }
        add = fault->block_size[0] * fault->block_size[1] * z + fault->block_size[0] * y + x;
        fault->param->block->ux2[add] = (float)uxm;
        fault->param->block->uy2[add] = (float)uym;
        fault->param->block->uz2[add] = (float)uzm;
        fault->param->block->crit2[add] = (float)critm;
        */

      }
}


static void fault_estimation_scale0_run(void *_fault)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  FAULT_ESTIMATION_PARAM *param = NULL;
  long x, y, z, n, offset;
  double c;

  if ( fault == NULL ) return;
  param = fault->param;
  if ( param == NULL ) return;
  offset = fault->param->block->border2;

  for ( z = 0; z < fault->param->block->scale0_size[2]; z++ )
    for ( y = 0; y < param->block->scale0_size[1]; y++ )
      for ( x = 0; x < param->block->scale0_size[0]; x++ )
        for ( n = 0; n < fault->array_phi0_size*fault->array_theta0_size; n++ )
        {
          c = pca_0(fault->crit_constant+1., param->block->gx, param->block->gy, param->block->gz, param->block->grad_size[0], param->block->grad_size[1], param->block->grad_size[2], x + offset, y + offset, z + offset,
                    param->X0[n], param->Y0[n], param->Z0[n], param->Npoint_per_disk0);
          param->block->scale0_crit[n][param->block->scale0_size[0] * param->block->scale0_size[1] * z + param->block->scale0_size[0] * y + x] = c;
        }
}

void fault_estimation_run(void *_fault, short *Ids)
{
  FAULT_ESTIMATION *fault = (FAULT_ESTIMATION*)_fault;
  FAULT_ESTIMATION_PARAM *param = NULL;
  FAULT_ESTIMATION_BLOCK_PARAM *block = NULL;
  long size, x, y, z, xx, yy, zz, b;
  float *tmp_f = NULL;

  if ( fault == NULL ) return;
  if ( fault->param == NULL )
    fault->param = fault_estimation_param3d_init(fault);
  if ( Ids == NULL ) return;
  param = fault->param;
  block = param->block;
  if ( block == NULL ) return;

	b = param->block->border1 + param->block->border2 + param->block->border3;
  size = block->block_size0[0] * block->block_size0[1] * block->block_size0[2];
	tmp_f = (float*)calloc(size, sizeof(float));
	
  for ( z = 0; z < block->block_size0[2]; z++ )
  {
    if ( z < b ) zz = b - z - 1; else if ( z > fault->depth - 1 + b ) zz = 2 * fault->depth - 1 + b - z; else zz = z - b;
    for ( y = 0; y < block->block_size0[1]; y++ )
    {
      if ( y < b ) yy = b - y - 1; else if ( y > fault->height - 1 + b ) yy = 2 * fault->height - 1 + b - y; else yy = y - b;
      for ( x = 0; x < block->block_size0[0]; x++ )
      {
        if ( x < b ) xx = b - x - 1; else if ( x > fault->width - 1 + b ) xx = 2 * fault->width - 1 + b - x; else xx = x - b;
        tmp_f[block->block_size0[0] * block->block_size0[1] * z + block->block_size0[0] * y + x] = (float)Ids[fault->width * fault->height * zz + fault->width * yy + xx];
      }    
    }
  }

  gradient_separable_valid_3d_float(CONVOLUTION_SEPARABLE_3D_VALID, tmp_f, block->block_size0[0], block->block_size0[1], block->block_size0[2],
                                    block->mask_grad_lp, block->mask_grad_hp, block->mask_grad_size, block->tmp1, block->tmp2, block->gx, block->gy, block->gz);
  fault_estimation_scale0_run(fault);
  fault_estimation_block_scale1_run(fault);
	FREE(tmp_f)
}


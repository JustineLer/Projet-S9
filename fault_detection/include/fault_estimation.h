// #pragma once

// #pragma once

#ifndef __FAULT_ESTIMATION__
#define __FAULT_ESTIMATION__

#include <config.h>

#define FAULT_ESTIMATION_CHRONOS_ALL            0
#define FAULT_ESTIMATION_CHRONOS_POSITION_CRIT  1
#define FAULT_ESTIMATION_CHRONOS_POSITION_CRIT0 2
#define FAULT_ESTIMATION_CHRONOS_TENSOR         3
#define FAULT_ESTIMATION_CHRONOS_TENSOR_CRIT    4
#define FAULT_ESTIMATION_CHRONOS_TENSOR_1       5
#define FAULT_ESTIMATION_CHRONOS_TENSOR_2       6
#define FAULT_ESTIMATION_CHRONOS_SCALE0         7
#define FAULT_ESTIMATION_CHRONOS_SCALE1         8
#define FAULT_ESTIMATION_CHRONOS_GRADIENT       9
#define FAULT_ESTIMATION_CHRONOS_WRITE_FILE    10
#define FAULT_ESTIMATION_CHRONOS_READ_FILE     11


#ifdef __cplusplus
extern "C" {
#endif

  EXPORT_LIB void *fault_estimation_init();

  EXPORT_LIB void *fault_estimation_release(void *_fault);

  EXPORT_LIB void fault_estimation_set_nbthreads(void *_fault, long nbthreads);

  EXPORT_LIB void fault_estimation_set_dims(void *_fault, long width, long height, long depth);
  EXPORT_LIB void fault_estimation_get_dims(void *_fault, long *width, long *height, long *depth);

  EXPORT_LIB void fault_estimation_set_array_rho0(void *_fault, double *rho, long size);
  EXPORT_LIB void fault_estimation_set_array_theta0(void *_fault, double *theta, long size);
  EXPORT_LIB void fault_estimation_set_array_phi0(void *_fault, double *phi, long size);
  EXPORT_LIB void fault_estimation_set_array_disk0(void *_fault, double *array_disk, long size);

  EXPORT_LIB void fault_estimation_set_array_rho1(void *_fault, double *theta, long size);
  EXPORT_LIB void fault_estimation_set_array_theta1(void *_fault, double *theta, long size);
  EXPORT_LIB void fault_estimation_set_array_phi1(void *_fault, double *phi, long size);
  EXPORT_LIB void fault_estimation_set_array_search1(void *_fault, double *pos, long size);

  EXPORT_LIB void fault_estimation_set_domain_in(void *_fault, double *X, double *Y, double *Z, long size);
  EXPORT_LIB void fault_estimation_get_domain_in(void *_fault, double **X, double **Y, double **Z, long *size);

  EXPORT_LIB void fault_estimation_set_out(void *_fault, float *nx, float *ny, float *nz, float *crit);

  EXPORT_LIB void fault_estimation_set_gradient(void *_fault, float *gx, float *gy, float *gz);
  // EXPORT_LIB void fault_estimation2_get_gradient(void *_fault, float **gx, float **gy, float **gz);

  EXPORT_LIB void fault_estimation_set_sigma_grad(void *_fault, double sigma);

  EXPORT_LIB void fault_estimation_set_data_in_filename(void *_fault, char *filename);
  EXPORT_LIB void fault_estimation_set_out_filename(void *_fault, char *nx_filename, char *ny_filename, char *nz_filename, char *crit_filename);
  EXPORT_LIB void fault_estimation_set_block_size(void *_fault, long bx, long by, long bz);

#define FAULT_ESTIMATION_RUN_TYPE_ALL   0
#define FAULT_ESTIMATION_RUN_TYPE_BLOCK 1
  EXPORT_LIB void fault_estimation_set_run_type(void *_fault, long type);
  EXPORT_LIB void fault_estimation_run(void *_fault, short *Ids);



#define FAULT_ESTIMATION_GET_DATA_OUT_NX    0
#define FAULT_ESTIMATION_GET_DATA_OUT_NY    1
#define FAULT_ESTIMATION_GET_DATA_OUT_NZ    2
#define FAULT_ESTIMATION_GET_DATA_OUT_CRIT  3
#define FAULT_ESTIMATION_GET_DATA_OUT_NX2   4
#define FAULT_ESTIMATION_GET_DATA_OUT_NY2   5
#define FAULT_ESTIMATION_GET_DATA_OUT_NZ2   6
#define FAULT_ESTIMATION_GET_DATA_OUT_CRIT2 7
  EXPORT_LIB float *fault_estimation_get_data_out(void *_fault, long type, long *width, long *height, long *depth);

  EXPORT_LIB void fault_estimation_set_chronos(void *_fault, long val);
  EXPORT_LIB void fault_estimation_chronos_reset(void *_fault, long label);
  EXPORT_LIB double fault_estimation_chronos_get(void *_fault, long label);

  // EXPORT_LIB void fault_estimation2_set_chronos(void *_fault, long val);

  // EXPORT_LIB void fault2_chronos_reset(void *_fault, long label);

  // EXPORT_LIB double fault2_estimation_chronos_get(void *_fault, long label);

  EXPORT_LIB void fault_estimation_set_crit_constant(void *_fault, double val);

  // EXPORT_LIB void fault_estimation2_set_normal_vector(void *_fault, double *v);

#define FAULT_ESTIMATION_MIX_TYPE_NONE   0
#define FAULT_ESTIMATION_MIX_TYPE_RAND   1
#define FAULT_ESTIMATION_MIX_TYPE_VECTOR 2
  // EXPORT_LIB void fault_estimation2_set_mix_type(void *_fault, long val);

  EXPORT_LIB long fault_estimation_get_conv_border(void *_fault);
  EXPORT_LIB void fault_estimation_get_nread_block(void *_fault, int *dim);

#ifdef __cplusplus
}
#endif

#endif
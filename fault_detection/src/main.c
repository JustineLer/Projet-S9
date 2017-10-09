
#include <stdio.h>
#include <malloc.h>

#include <util.h>
#include <signal.h>
#include <fileIO.h>
#include <chronos.h>
#include <fault_estimation.h>

int main(int argc, char **argv)
{

  void *fault = NULL, *chronos = NULL;
  double sigma = 1.0, *X = NULL, *Y = NULL, *Z = NULL, *array_theta0 = NULL, *array_phi0 = NULL, *array_theta1 = NULL, *array_phi1 = NULL, *array_rho0 = NULL, *array_rho1 = NULL,
    *array_disk = NULL, *array_pos = NULL, *array_disk0 = NULL, *array_search1 = NULL,
    rho0_min = 0.0, rho0_max = 3.0,
    theta0_min = 0.0, theta0_max = 2.0*PI,
    phi0_min = 0.0, phi0_max = PI,
    array_disk0_min = -1.0, array_disk0_max = 1.0,
    rho1_min = 0.0, rho1_max = 9.0,
    theta1_min = 0.0, theta1_max = 2.0*PI,
    phi1_min = PI / 4, phi1_max = PI / 2.0,
    array_search1_min = -2.0, array_search1_max = 2.0;
  float *h_hp = NULL, *h_lp = NULL, *temp1 = NULL, *temp2 = NULL, *gx = NULL, *gy = NULL, *gz = NULL, *data = NULL, *nx = NULL, *ny = NULL, *nz = NULL, *crit = NULL,
    *nx_out = NULL, *ny_out = NULL, *nz_out = NULL, *crit_out = NULL;
  long size,
    Nrho0 = 4, Ntheta0 = 8, Nphi0 = 4, array_disk0_size = 3, Nrho1 = 3, Ntheta1 = 30, Nphi1 = 4, array_search1_size = 5, Nthreads = 6;
  int width, height, depth;

  char in_filename[1000] = "data128_int16.gdr",
    nx_filename[1000] = "nx0.gdr",
    ny_filename[1000] = "ny0.gdr",
    nz_filename[1000] = "nz0.gdr",
    crit_filename[1000] = "crit0.gdr";
  double t_total/*, tscale0, tscale1, tgradient, twrite, tread*/;
  short *Ids = NULL, *Ids0 = NULL;

  chronos = CHRONOS_INIT;

  fileGetDims(in_filename, &width, &height, &depth);
  size = width*height*depth;

  array_rho0 = (double*)calloc(Nrho0, sizeof(double));
  linspace_fill(rho0_min, rho0_max, Nrho0, array_rho0);

  array_theta0 = (double*)calloc(Ntheta0, sizeof(double));
  linspace_fill(theta0_min, theta0_max*(1.0 - 1.0 / Ntheta0), Ntheta0, array_theta0);

  array_phi0 = (double*)calloc(Nphi0, sizeof(double));
  linspace_fill(phi0_min, phi0_max*(1.0 - 1.0 / Nphi0), Nphi0, array_phi0);

  array_disk0 = (double*)calloc(array_disk0_size, sizeof(double));
  linspace_fill(array_disk0_min, array_disk0_max, array_disk0_size, array_disk0);

  array_rho1 = (double*)calloc(Nrho1, sizeof(double));
  linspace_fill(rho1_min, rho1_max, Nrho1, array_rho1);

  array_theta1 = (double*)calloc(Ntheta1, sizeof(double));
  linspace_fill(theta1_min, theta1_max*(1.0 - 1.0 / Ntheta1), Ntheta1, array_theta1);

  array_phi1 = (double*)calloc(Nphi1, sizeof(double));
  linspace_fill(phi1_min, phi1_max, Nphi1, array_phi1);

  array_search1 = (double*)calloc(array_search1_size, sizeof(double));
  linspace_fill(array_search1_min, array_search1_max, array_search1_size, array_search1);

  fault = fault_estimation_init();
  fault_estimation_set_nbthreads(fault, Nthreads);
  fault_estimation_set_block_size(fault, width, height, depth);
  fault_estimation_set_dims(fault, width, height, depth);

  fault_estimation_set_array_theta0(fault, array_theta0, Ntheta0);
  fault_estimation_set_array_phi0(fault, array_phi0, Nphi0);
  fault_estimation_set_array_rho0(fault, array_rho0, Nrho0);
  fault_estimation_set_array_disk0(fault, array_disk0, array_disk0_size);

  fault_estimation_set_array_theta1(fault, array_theta1, Ntheta1);
  fault_estimation_set_array_phi1(fault, array_phi1, Nphi1);
  fault_estimation_set_array_rho1(fault, array_rho1, Nrho1);
  fault_estimation_set_array_search1(fault, array_search1, array_search1_size);

  Ids0 = (short*)fileReadBlockByName(in_filename);
  CHRONOS_TIC(chronos, 0);
  fault_estimation_run(fault, Ids0);
  CHRONOS_TOC(chronos, 0);
  t_total = CHRONOS_GET(chronos, 0);  
	Ids0 = fileFreeData(Ids0);
  
  nx_out = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_NX, NULL, NULL, NULL);
  ny_out = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_NY, NULL, NULL, NULL);
  nz_out = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_NZ, NULL, NULL, NULL);
  crit_out = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_CRIT, NULL, NULL, NULL);

  fileWriteBlock(nx_filename, GDR_FORMAT_FLOAT, nx_out, width, height, depth);
  fileWriteBlock(ny_filename, GDR_FORMAT_FLOAT, ny_out, width, height, depth);
  fileWriteBlock(nz_filename, GDR_FORMAT_FLOAT, nz_out, width, height, depth);
  fileWriteBlock(crit_filename, GDR_FORMAT_FLOAT, crit_out, width, height, depth);

  FREE(array_rho0)
  FREE(array_theta0)
  FREE(array_phi0)
  FREE(array_disk0)
  FREE(array_rho1)
  FREE(array_theta1)
  FREE(array_phi1)
  FREE(array_search1)
  fault = fault_estimation_release(fault);

  printf("temps total: %.1f ms\n", t_total);
  printf("\nend.\n");
  return 1;
}

#ifndef _FAULT_ESTIMATION_STRUCT_
#define _FAULT_ESTIMATION_STRUCT_


typedef struct _FAULT_ESTIMATION_BLOCK_PARAM
{
  long block_size0[3], grad_size[3], scale0_size[3], border1, border2, border3, mask_grad_size;
  float *ux, *uy, *uz, *crit, *mask_grad_lp, *mask_grad_hp, *ux2, *uy2, *uz2, *crit2;
  double **scale0_crit;
} FAULT_ESTIMATION_BLOCK_PARAM;


typedef struct _FAULT_ESTIMATION_PARAM
{
  long ** X0, ** Y0, ** Z0;
  long ** X1, ** Y1, ** Z1;
  long Nangle0, Npoint_per_disk0;
  long Nangle1, Npoint_per_disk1;
  long plane_border[3],   // plus grand pour la recherche des plans
       grad_border,       // plus grand pour le calcul du gradient
       border[3],         // bordure totale (plan + gradient)
       * scale1_to_scale0_angle_match;
  double * ux0, * uy0, * uz0;
	double * ux1, * uy1, * uz1;
  int ** array_exclusion;
  FAULT_ESTIMATION_BLOCK_PARAM block;
} FAULT_ESTIMATION_PARAM;


typedef struct _FAULT_ESTIMATION
{
  long width, height, depth;
  long array_rho0_size, array_theta0_size, array_phi0_size, array_ortho0_size;
  long array_rho1_size, array_theta1_size, array_phi1_size, array_ortho1_size;
  long size_in, block_size[3];
  double * array_rho0, * array_theta0, * array_phi0, * array_ortho0;
  double * array_rho1, * array_theta1, * array_phi1, * array_ortho1; 
  double crit_constant, sigma_grad;
	float * gx, * gy, * gz;  
	float * nx, * ny, * nz;
  float * crit;
  FAULT_ESTIMATION_PARAM * param;
} FAULT_ESTIMATION;

#endif
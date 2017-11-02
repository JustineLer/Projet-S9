
#include <stdlib.h>

#include <fault_estimation_interface.h>
#include <fault_estimation_struct.h>

#include <fault_estimation.h>
#include <util.h>


void * fault_estimation_init()
{
  FAULT_ESTIMATION * fault_estimation = NULL;

  fault_estimation = (FAULT_ESTIMATION *)calloc(1, sizeof(FAULT_ESTIMATION));
  fault_estimation_set_block_size(fault_estimation, 64, 64, 64);
  fault_estimation_set_sigma_grad(fault_estimation, 1.0);
  return (void *)fault_estimation;
}


static void fault_estimate_block_param_realease(FAULT_ESTIMATION_BLOCK_PARAM * block)
{
  if ( block == NULL ) return;
  FREE_VTAB(block->scale0_crit)
  FREE(block->ux)
  FREE(block->uy)
  FREE(block->uz)
  FREE(block->crit)
  FREE(block->mask_grad_lp)
  FREE(block->mask_grad_hp)
  FREE_VTAB(block->scale0_crit)
}


static FAULT_ESTIMATION_PARAM * fault_estimation_param_release(FAULT_ESTIMATION_PARAM * param)
{
  if ( param == NULL ) return NULL;
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
  fault_estimate_block_param_realease(&(param->block));
  FREE(param)
  return NULL;
}


void * fault_estimation_release(void * fault_estimation)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;

  if ( fault_estimation == NULL ) return NULL;
  _fault_estimation->param = fault_estimation_param_release(_fault_estimation->param);
  FREE(_fault_estimation->array_rho0)
  FREE(_fault_estimation->array_theta0)
  FREE(_fault_estimation->array_phi0)
	FREE(_fault_estimation->array_ortho0)
	FREE(_fault_estimation->array_rho1)
  FREE(_fault_estimation->array_theta1)
  FREE(_fault_estimation->array_phi1)
  FREE(_fault_estimation->array_ortho1)
  FREE(_fault_estimation)
  return NULL;
}


static void linspace(double x1, double x2, long nb, double * v)
{
  long n;
  
  for (n=0; n<nb; n++) v[n] = x1+(x2-x1)/(double)(nb-1)*(double)n;
}


static void set_linspace(double ** ptr_data, long * ptr_nb, double min, double max, long nb)
{
	FREE(*ptr_data)
  *ptr_data = (double *)calloc(nb, sizeof(double));
  *ptr_nb = nb;
	linspace(min, max, nb, *ptr_data);
}


void fault_estimation_set_rho0(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_rho0, &_fault_estimation->array_rho0_size, min, max, nb);
}


void fault_estimation_set_theta0(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_theta0, &_fault_estimation->array_theta0_size, min, max, nb);
}


void fault_estimation_set_phi0(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_phi0, &_fault_estimation->array_phi0_size, min, max, nb);
}


void fault_estimation_set_ortho0(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_ortho0, &_fault_estimation->array_ortho0_size, min, max, nb);
}


void fault_estimation_set_rho1(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_rho1, &_fault_estimation->array_rho1_size, min, max, nb);
}


void fault_estimation_set_theta1(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_theta1, &_fault_estimation->array_theta1_size, min, max, nb);
}


void fault_estimation_set_phi1(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_phi1, &_fault_estimation->array_phi1_size, min, max, nb);
}


void fault_estimation_set_ortho1(void * fault_estimation, double min, double max, long nb)
{
  FAULT_ESTIMATION * _fault_estimation = (FAULT_ESTIMATION *)fault_estimation;
 
  if ( fault_estimation != NULL )
		set_linspace(&_fault_estimation->array_ortho1, &_fault_estimation->array_ortho1_size, min, max, nb);
}

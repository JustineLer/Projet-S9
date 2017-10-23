
#include <malloc.h>
#include <string.h>
#include <math.h>

#include <callback.h>
#include <util.h>
#include <signal.h>
#include <chronos.h>
#include <fault_estimation.h>


#define PROG_NAME "Fault orientatoion"
#define NUM_VERSION  "0.1"
#define DATE_VERSION __DATE__

static void On_Menu_Variable_Size()
{
	// ndMsgFormat("size bool %d", 8*sizeof(bool));
	ndMsgFormat("size char %d", 8*sizeof(char)); 
	ndMsgFormat("size short %d", 8*sizeof(short)); 
	ndMsgFormat("size int %d", 8*sizeof(int)); 
	ndMsgFormat("size long %d", 8*sizeof(long)); 
	ndMsgFormat("size long long %d", 8*sizeof(long long)); 
	ndMsgFormat("size float %d", 8*sizeof(float)); 
	ndMsgFormat("size double %d", 8*sizeof(double)); 
	ndMsgFormat("size long double %d", 8*sizeof(long double));
	ndMsgFormat("size __int8 %d", 8*sizeof(__int8));
	ndMsgFormat("size __int16 %d", 8*sizeof(__int16));
	ndMsgFormat("size __int32 %d", 8*sizeof(__int32));
	ndMsgFormat("size __int64 %d", 8*sizeof(__int64));
	ndMsgFormat("size wchar_t %d", 8*sizeof(wchar_t));
	ndMsgFormat("size size_t %d", 8*sizeof(size_t));
	ndMsgFormat("size pointeur %d", 8*sizeof(void*));
}

static void On_Menu_Version()
{
	ndMsgBox(PROG_NAME, "version: "NUM_VERSION"\n\ndate: "DATE_VERSION);
}


static void On_Menu_fault_orientation_cpu_3D()
{
  long id_in = 0, id_nx = 0, id_ny = 0, id_nz = 0, id_crit = 0,    
    width, height, depth, Nrho0 = 4, Ntheta0 = 8, Nphi0 = 4, array_disk0_size = 3, Nrho1 = 3, Ntheta1 = 30, Nphi1 = 4, array_search1_size = 5, Nthreads = 4, z, add, n, id;
  double rho0_min = 0.0, rho0_max = 3.0,
    theta0_min = 0.0, theta0_max = 2.0*PI,
    phi0_min = 0.0, phi0_max = PI,
    array_disk0_min = -1.0, array_disk0_max = 1.0,
    rho1_min = 0.0, rho1_max = 9.0,
    theta1_min = 0.0, theta1_max = 2.0*PI,
    phi1_min = PI / 4, phi1_max = PI / 2.0,
    array_search1_min = -2.0, array_search1_max = 2.0,
    *array_rho0 = NULL, *array_theta0 = NULL, *array_phi0 = NULL, *array_disk0 = NULL, *array_rho1 = NULL, *array_theta1 = NULL, *array_phi1 = NULL, *array_search1 = NULL,
    t_total = 0.0;
  void * box = NULL, *fault = NULL, *chronos = NULL;
  short *Ids = NULL, *tmp = NULL;
  float *ptr = NULL;

  box = ndBoxParameterAdd(box, "scale 0%t%n");
  box = ndBoxParameterAdd(box, "rho0 min%t%f%trho0 max%f%tNrho0%d%n", &rho0_min, &rho0_max, &Nrho0);
  box = ndBoxParameterAdd(box, "theta0 min%t%f%ttheta0 max%f%tNtheta0%d%n", &theta0_min, &theta0_max, &Ntheta0);
  box = ndBoxParameterAdd(box, "phi0 min%t%f%tphi0 max%f%tNphi0%d%n", &phi0_min, &phi0_max, &Nphi0);
  box = ndBoxParameterAdd(box, "array disk0 min%t%f%tarray disk0 max%f%tNarray disk0%d%n", &array_disk0_min, &array_disk0_max, &array_disk0_size);

  box = ndBoxParameterAdd(box, "scale 1%t%n");
  box = ndBoxParameterAdd(box, "rho1 min%f%trho1 max%f%tNrho1%d%n", &rho1_min, &rho1_max, &Nrho1);
  box = ndBoxParameterAdd(box, "theta1 min%f%ttheta1 max%f%tNtheta1%d%n", &theta1_min, &theta1_max, &Ntheta1);
  box = ndBoxParameterAdd(box, "phi1 min%f%tphi1 max%f%tNphi1%d%n", &phi1_min, &phi1_max, &Nphi1);
  box = ndBoxParameterAdd(box, "array search min%f%tarray search max%f%tN array search%d%n", &array_search1_min, &array_search1_max, &array_search1_size);
  
  if ( !ndBoxParameterDialog(box) ) return;
  // =======================================
  id_in = ndBlockGetIn();
  ndBlockGetSize(id_in, &width, &height, &depth);
  Ids = (short*)calloc(width*height*depth, sizeof(short));
  for ( z = 0; z < depth; z++ )
  {
    tmp = ndBlockGetFront(id_in, z, TYPE_SHORT);
    for ( add = 0; add < width*height; add++ )
      Ids[width*height*z + add] = tmp[add];
    ndFree(tmp);
  }

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
  // scale 0
  fault_estimation_set_array_theta0(fault, array_theta0, Ntheta0);
  fault_estimation_set_array_phi0(fault, array_phi0, Nphi0);
  fault_estimation_set_array_rho0(fault, array_rho0, Nrho0);
  fault_estimation_set_array_disk0(fault, array_disk0, array_disk0_size);
  // scale 1
  fault_estimation_set_array_theta1(fault, array_theta1, Ntheta1);
  fault_estimation_set_array_phi1(fault, array_phi1, Nphi1);
  fault_estimation_set_array_rho1(fault, array_rho1, Nrho1);
  fault_estimation_set_array_search1(fault, array_search1, array_search1_size);

  ndMsg("CPU RUN");
  chronos = CHRONOS_INIT;
  CHRONOS_TIC(chronos, 0);
  fault_estimation_run(fault, Ids);
  CHRONOS_TOC(chronos, 0);
  t_total = CHRONOS_GET(chronos, 0);
  FREE(Ids)

  id_nx = ndBlockNew();
  id_ny = ndBlockNew();
  id_nz = ndBlockNew();
  id_crit = ndBlockNew();
  for ( n = 0; n < 4; n++ )
  {
    switch ( n )
    {
      case 0:
        id = id_nx; 
        ptr = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_NX, NULL, NULL, NULL);
        break;
      case 1:
        id = id_ny;
        ptr = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_NY, NULL, NULL, NULL);
        break;
      case 2:
        id = id_nz;
        ptr = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_NZ, NULL, NULL, NULL);
        break;
      case 3:
        id = id_crit;
        ptr = fault_estimation_get_data_out(fault, FAULT_ESTIMATION_GET_DATA_OUT_CRIT, NULL, NULL, NULL);
        break;
    }
    ndBlockSetTypeData(id, TYPE_FLOAT);
    ndBlockSetSize(id, width, height, depth);
    ndBlockShow(id);
    if ( ptr != NULL )
      for ( z = 0; z<depth; z++ )
        ndBlockSetFront(id, z, ptr + z*width*height);
  }
  FREE(array_rho0)
  FREE(array_theta0)
  FREE(array_phi0)
  FREE(array_disk0)
  FREE(array_rho1)
  FREE(array_theta1)
  FREE(array_phi1)
  FREE(array_search1)
  fault = fault_estimation_release(fault);
  ndMsgFormat("cpu fault estimation total time: %.1f ms", t_total);
}

static void On_Menu_fault_orientation_gpu_3d()
{

}

EXPORT_LIB void attLibraryInit(long id_library)
{
	long menu;

  menu = ndMenuCreate("Fault Orientation");
  ndMenuAddItem(menu, "Fault orientation (CPU)", On_Menu_fault_orientation_cpu_3D);
  ndMenuAddItem(menu, "Fault Orientation (GPU)", On_Menu_fault_orientation_gpu_3d);
	ndMenuAddItem(menu, "Size Variables", On_Menu_Variable_Size);
	ndMenuAddItem(menu, "Version", On_Menu_Version);
  ndLibraryAddMenu(id_library, "MENU_3D", menu);	
}


EXPORT_LIB void attLibraryEnd(long id_library)
{

}

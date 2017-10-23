
#include <malloc.h>
#include <string.h>
#include <math.h>

#include <callback.h>
#include <util.h>
#include <signal.h>
#include <chronos.h>
#include <fault_estimation_interface.h>
#include <fault_estimation.h>


#define PROG_NAME "Fault"
#define NUM_VERSION  "0.2"
#define DATE_VERSION __DATE__


double rho0_min = 0.0, rho0_max = 3.0,
			 theta0_min = 0.0, theta0_max = 360.0,
       phi0_min = 0.0, phi0_max = 180.0,
			 ortho0_min = -1.0, ortho0_max = 1.0,
       rho1_min = 0.0, rho1_max = 9.0,
       theta1_min = 0.0, theta1_max = 360.0,
       phi1_min = 45.0, phi1_max = 90.0,
			 ortho1_min = -2.0, ortho1_max = 2.0;
long rho0_nb = 4, theta0_nb = 8, phi0_nb = 4, ortho0_nb = 3, rho1_nb = 3, theta1_nb = 30, phi1_nb = 4, ortho1_nb = 5;


static void On_Menu_Fault_Parameters()
{
	void * box = NULL;
	
	box = ndBoxParameterAdd(box, "Echelle 0%t%n");
  box = ndBoxParameterAdd(box, "Rayon%tMin%fMax%fNb%d%n", &rho0_min, &rho0_max, &rho0_nb);
  box = ndBoxParameterAdd(box, "Theta(°)%tMin%fMax%fNb%d%n", &theta0_min, &theta0_max, &theta0_nb);
  box = ndBoxParameterAdd(box, "Phi(°)%tMin%fMax%fNb%d%n", &phi0_min, &phi0_max, &phi0_nb);
  box = ndBoxParameterAdd(box, "Orthogonal%tMin%fMax%fNb%d%n", &ortho0_min, &ortho0_max, &ortho0_nb);
  box = ndBoxParameterAdd(box, "Echelle 1%t%n");
  box = ndBoxParameterAdd(box, "Rayon%tMin%fMax%fNb%d%n", &rho1_min, &rho1_max, &rho1_nb);
  box = ndBoxParameterAdd(box, "Theta(°)%tMin%fMax%fNb%d%n", &theta1_min, &theta1_max, &theta1_nb);
  box = ndBoxParameterAdd(box, "Phi(°)%tMin%fMax%fNb%d%n", &phi1_min, &phi1_max, &phi1_nb);
	box = ndBoxParameterAdd(box, "Orthogonal%tMin%fMax%fNb%d%n", &ortho1_min, &ortho1_max, &ortho1_nb);
  ndBoxParameterDialog(box);
}


static float * in_to_data_float(long id, long x1, long x2, long y1, long y2, long z1, long z2)
{
	long w, h, d, z, dx = x2+1-x1, dy = y2+1-y1, dz = z2+1-z1, add;
	float * data = NULL, * tmp = NULL;

	ndBlockGetSize(id, &w, &h, &d);
  data = (float *)calloc(dx*dy*dz, sizeof(float));
  for (z=z1; z<=z2; z++)
		{
			tmp = (float *)ndBlockGetPartialFront(id, z, TYPE_FLOAT, x1, y1, dx, dy);
			for (add=0; add<dx*dy; add++) data[dx*dy*(z-z1)+add] = tmp[add];
			ndFree(tmp);
		}
	return data;
}


static void data_float_to_out(float * ptr, long width, long height, long depth, char * name)
{
	long id, z;
 
	id = ndBlockNew();
	ndBlockSetTypeData(id, TYPE_FLOAT);
  ndBlockSetSize(id, width, height, depth);
  ndBlockShow(id);
  if ( ptr != NULL )
		for (z=0; z<depth; z++) ndBlockSetFront(id, z, ptr+z*width*height);
	ndWindowSetName(id, name);
}


static void On_Menu_Fault_CPU_GPU(long type)
{
  long id_in, width, height, depth;	
  double dt = 0.0;
  void * fault_estimation = NULL, * chronos = NULL;
  float * data = NULL;

  id_in = ndBlockGetIn();
  ndBlockGetSize(id_in, &width, &height, &depth);
	data = in_to_data_float(id_in, 0, width-1, 0, height-1, 0, depth-1);

  fault_estimation = fault_estimation_init();
  fault_estimation_set_block_size(fault_estimation, width, height, depth);
  fault_estimation_set_dims(fault_estimation, width, height, depth);
  fault_estimation_set_theta0(fault_estimation, theta0_min, theta0_max, theta0_nb);
  fault_estimation_set_phi0(fault_estimation, phi0_min, phi0_max, phi0_nb);
  fault_estimation_set_rho0(fault_estimation, rho0_min, rho0_max, rho0_nb);
  fault_estimation_set_ortho0(fault_estimation, ortho0_min, ortho0_max, ortho0_nb);
  fault_estimation_set_theta1(fault_estimation, theta1_min, theta1_max, theta1_nb);
  fault_estimation_set_phi1(fault_estimation, phi1_min, phi1_max, phi1_nb);
  fault_estimation_set_rho1(fault_estimation, rho1_min, rho1_max, rho1_nb);
  fault_estimation_set_ortho1(fault_estimation, ortho1_min, ortho1_max, ortho1_nb); 

	if ( type == 0 ) ndMsg("CPU RUN"); else ndMsg("GPU RUN");

  chronos = CHRONOS_INIT;
  CHRONOS_TIC(chronos, 0);
  fault_estimation_run(fault_estimation, data);
  CHRONOS_TOC(chronos, 0);
  dt = CHRONOS_GET(chronos, 0);
	ndMsgFormat("Temps d'exécution : %.1f ms", dt);

  free(data);
	data = fault_estimation_get_data_out(fault_estimation, FAULT_ESTIMATION_GET_DATA_OUT_NX, NULL, NULL, NULL);
	data_float_to_out(data, width, height, depth, "Nx");
	data = fault_estimation_get_data_out(fault_estimation, FAULT_ESTIMATION_GET_DATA_OUT_NY, NULL, NULL, NULL);
	data_float_to_out(data, width, height, depth, "Ny");
	data = fault_estimation_get_data_out(fault_estimation, FAULT_ESTIMATION_GET_DATA_OUT_NZ, NULL, NULL, NULL);
	data_float_to_out(data, width, height, depth, "Nz");
	data = fault_estimation_get_data_out(fault_estimation, FAULT_ESTIMATION_GET_DATA_OUT_CRIT, NULL, NULL, NULL);
	data_float_to_out(data, width, height, depth, "Critère");
  fault_estimation = fault_estimation_release(fault_estimation);
}


static void On_Menu_Fault_CPU()
{
	On_Menu_Fault_CPU_GPU(0);
}


static void On_Menu_Fault_GPU()
{
	On_Menu_Fault_CPU_GPU(1);
}


static void On_Menu_Fault_Version()
{
	ndMsgBox(PROG_NAME, "version: "NUM_VERSION"\n\ndate: "DATE_VERSION);
}


EXPORT_LIB void attLibraryInit(long id_library)
{
	long menu;

  menu = ndMenuCreate("Détection de failles");
	ndMenuAddItem(menu, "Paramètres", On_Menu_Fault_Parameters);
	ndMenuAddSeparator(menu);
  ndMenuAddItem(menu, "Version CPU", On_Menu_Fault_CPU);
  ndMenuAddItem(menu, "Version GPU", On_Menu_Fault_GPU);
	ndMenuAddSeparator(menu);
	ndMenuAddItem(menu, "Version", On_Menu_Fault_Version);
  ndLibraryAddMenu(id_library, "MENU_3D", menu);	
}


EXPORT_LIB void attLibraryEnd(long id_library)
{

}

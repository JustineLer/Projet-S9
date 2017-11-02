
void * fault_estimation_init();
void * fault_estimation_release(void * fault_estimation);

void fault_estimation_set_rho0(void * fault_estimation, double min, double max, long nb);
void fault_estimation_set_theta0(void * fault_estimation, double min, double max, long nb);
void fault_estimation_set_phi0(void * fault_estimation, double min, double max, long nb);
void fault_estimation_set_ortho0(void * fault_estimation, double min, double max, long nb);
void fault_estimation_set_rho1(void * fault_estimation, double min, double max, long nb);
void fault_estimation_set_theta1(void * fault_estimation, double min, double max, long nb);
void fault_estimation_set_phi1(void * fault_estimation, double min, double max, long nb);
void fault_estimation_set_ortho1(void * fault_estimation, double min, double max, long nb);

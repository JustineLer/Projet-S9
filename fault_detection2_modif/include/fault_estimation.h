#ifndef __FAULT_ESTIMATION__
#define __FAULT_ESTIMATION__

#ifdef __cplusplus
extern "C" {
#endif

void fault_estimation_set_dims(void * fault_estimation, long width, long height, long depth);
void fault_estimation_get_dims(void * fault_estimation, long * width, long * height, long * depth);

void fault_estimation_set_sigma_grad(void * fault_estimation, double sigma);

void fault_estimation_set_block_size(void * fault_estimation, long bx, long by, long bz);

void fault_estimation_run(void * fault_estimation, float * data);

#define FAULT_ESTIMATION_GET_DATA_OUT_NX    0
#define FAULT_ESTIMATION_GET_DATA_OUT_NY    1
#define FAULT_ESTIMATION_GET_DATA_OUT_NZ    2
#define FAULT_ESTIMATION_GET_DATA_OUT_CRIT  3
float * fault_estimation_get_data_out(void * fault_estimation, long type, long * width, long * height, long * depth);

#ifdef __cplusplus
}
#endif

#endif
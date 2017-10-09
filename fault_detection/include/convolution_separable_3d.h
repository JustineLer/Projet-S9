
#ifndef __CONVOLUTION_SEPARABLE_3D__
#define __CONVOLUTION_SEPARABLE_3D__
#endif

#define CONVOLUTION_SEPARABLE_3D_VALID 0
#define CONVOLUTION_SEPARABLE_3D_SAME  1

#define CONVOLUTION_SEPARABLE_3D_DOUBLE 0
#define CONVOLUTION_SEPARABLE_3D_FLOAT  1


#ifdef __cplusplus
extern "C" {
#endif

	EXPORT_LIB void gradient_separable_valid_3d_float(long type, float *data, long width, long height, long depth,
    float *h_lp, float *h_hp, long h_size,
    float *temp1, float *temp2,
    float *gx, float *gy, float *gz);


#ifdef __cplusplus
}
#endif
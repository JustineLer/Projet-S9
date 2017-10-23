
#ifdef __cplusplus
extern "C" {
#endif
	void gradient_gpu_main(float * img_in, int width, int height, int depth,
		float *filterSmooth, float *derivate, int t, float *gx, float*gy, float *gz);

#ifdef __cplusplus
}
#endif

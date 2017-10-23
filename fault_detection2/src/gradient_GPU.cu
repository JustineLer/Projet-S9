#include<cuda.h>
#include<cuda_runtime.h>
#include"gradient_GPU.cuh"



__global__ void convX_gpu(float *img_in, int width, int height, int depth, float *filter, int t, float *outX)
{
	int X = blockIdx.x*blockDim.x + threadIdx.x;
	int Y = blockIdx.y*blockDim.y + threadIdx.y;
	int Z = blockIdx.z*blockDim.z + threadIdx.z;
	int indice;
	float s = 0.0;

	if (X > t / 2 && X < width - t / 2 && Y >= 0 && Y < height && Z >= 0 && Z < depth) 
	{
		for (int j = -t / 2; j <= t / 2; j++)
		{
			indice = width*(height*Z + Y) + j + X;
			s += img_in[indice] * filter[t / 2 - j];
		}

		outX[width*(height*Z + Y) + X] = (float)s;
	}
}


__global__ void convY_gpu(float *img_in, int width, int height, int depth, float *filter, int t, float *outY)
{
	int X = blockIdx.x*blockDim.x + threadIdx.x;
	int Y = blockIdx.y*blockDim.y + threadIdx.y;
	int Z = blockIdx.z*blockDim.z + threadIdx.z;
	int indice;
	float s = 0.0;

	if (Y > t / 2 && Y < height - t / 2 && X >= 0 && X < width && Z >= 0 && Z < depth)
	{
		for (int j = -t / 2; j <= t / 2; j++)
		{
			indice = width*(height*Z + j+Y) + X;
			s += img_in[indice] * filter[t / 2 - j];
		}

		outY[width*(height*Z + Y) + X] = (float)s;
	}
}
__global__ void convZ_gpu(float *img_in, int width, int height, int depth, float *filter, int t, float *outZ)
{
	int X = blockIdx.x*blockDim.x + threadIdx.x;
	int Y = blockIdx.y*blockDim.y + threadIdx.y;
	int Z = blockIdx.z*blockDim.z + threadIdx.z;
	int indice;
	float s = 0.0;

	if (Z > t / 2 && Z < depth - t / 2 && Y >= 0 && Y < height && X >= 0 && X < width)
	{
		for (int j = -t / 2; j <= t / 2; j++)
		{
			indice = width*(height*(Z+j) + Y) + X;
			s += img_in[indice] * filter[t / 2 - j];
		}

		outZ[width*(height*Z + Y) + X] = (float)s;
	}
}



void gradient_gpu_main(float * img_in, int width, int height, int depth,
	float *filterSmooth, float *derivate, int t, float *gx, float*gy, float *gz)
{

	int nbThreadX = 10, nbThreadY = 10, nbThreadZ = 10;
	int size = width*height*depth;
	dim3 dimBlock(nbThreadX,nbThreadY,nbThreadZ);
	dim3 dimGrid((width + nbThreadX - 1) / nbThreadX, (height + nbThreadY - 1) / nbThreadY, (depth + nbThreadZ - 1) / nbThreadZ);
	float *img_in_cuda = NULL;
	float *filterSmooth_cuda = NULL;
	float *derivate_cuda = NULL;
	float *gx_cuda = NULL, *gy_cuda = NULL, *gz_cuda = NULL;
	float *tmp1, *tmp2;

	cudaMalloc((void **)&img_in_cuda, size *sizeof(float));
	cudaMemcpy(img_in_cuda, img_in, size * sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void **)&filterSmooth_cuda, t* sizeof(float));
	cudaMemcpy(filterSmooth_cuda, filterSmooth, t * sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void **)&derivate_cuda, t * sizeof(float));
	cudaMemcpy(derivate_cuda, derivate, t * sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void **)&gx_cuda, size * sizeof(float));
	cudaMalloc((void **)&gy_cuda, size * sizeof(float));
	cudaMalloc((void **)&gz_cuda, size * sizeof(float));

	cudaMalloc((void **)&tmp1, size * sizeof(float));
	cudaMalloc((void **)&tmp2, size * sizeof(float));

	convZ_gpu <<<dimGrid, dimBlock >>> (img_in_cuda, width, height, depth, filterSmooth_cuda, t, tmp1);
	cudaThreadSynchronize();
	convY_gpu <<<dimGrid, dimBlock >>> (tmp1, width, height, depth, filterSmooth_cuda, t, tmp2);
	cudaThreadSynchronize();
	convX_gpu <<<dimGrid, dimBlock >>> (tmp2, width, height, depth, derivate_cuda, t, gx_cuda);
	cudaThreadSynchronize();
	
	convX_gpu <<<dimGrid, dimBlock >>> (tmp1, width, height, depth, filterSmooth_cuda, t, tmp2);
	cudaThreadSynchronize();
	convY_gpu <<<dimGrid, dimBlock >>> (tmp2, width, height, depth, derivate_cuda, t, gy_cuda);
	cudaThreadSynchronize();

	convX_gpu<<<dimGrid, dimBlock >>> (img_in_cuda, width, height, depth, filterSmooth_cuda, t, tmp1);
	cudaThreadSynchronize();
	convY_gpu<<<dimGrid, dimBlock >>> (tmp1, width, height, depth, filterSmooth_cuda, t, tmp2);
	cudaThreadSynchronize(); 
	convZ_gpu<<<dimGrid, dimBlock >>> (tmp2, width, height, depth, derivate_cuda, t, gz_cuda);
	cudaThreadSynchronize();

	
	cudaMemcpy(gx_cuda, gx, size * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(gy_cuda, gy, size * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(gz_cuda, gz, size * sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(gx_cuda);
	cudaFree(gy_cuda);
	cudaFree(gz_cuda);

	cudaFree(img_in_cuda);
	cudaFree(filterSmooth_cuda);
	cudaFree(derivate_cuda);
	cudaFree(tmp1);
	cudaFree(tmp2);

}

 
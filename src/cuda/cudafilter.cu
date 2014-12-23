/*
 * cudafilter.cu
 *
 *  Created on: 9 Jan, 2013
 *      Author: zhao
 */
#define MASK_WIDTH 3

#include <stdio.h>
#include <cuda.h>
#include "../util5/cudpp.h"
#include "../util5/helper_cuda.h"
#include "cudatypes.h"

__constant__ float c_mask[MASK_WIDTH*MASK_WIDTH];
/**
 * we assume that mask_width/2 is smaller than the block width
 * we assume that mask_width = mask_height
 * if it is not the case, the kernel should be reimplemented.
 */
__global__ void convolve_kernel(float* d_pixArrayIn,
		float* d_pixArrayOut,
		unsigned int width,
		unsigned int height,
		unsigned int mask_width) {

	const int id_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int id_y = blockDim.y * blockIdx.y + threadIdx.y;

	//the half width of the convolve mask(the margin width of the tiling area)
	int half_mask_width = mask_width / 2;

	//in case of the margin width is larger than the block width
	//int n = half_mask_width / SQUARE_BLK_WIDTH + 1;

	__shared__ float s_pix[(SQUARE_BLK_WIDTH + MASK_WIDTH - 1)][(SQUARE_BLK_WIDTH + MASK_WIDTH - 1)];

	//load the upper left area
	int halo_index_x = (blockIdx.x-1)*blockDim.x + threadIdx.x;
	int halo_index_y = (blockIdx.y-1)*blockDim.y + threadIdx.y;

	if(threadIdx.x >= blockDim.x - half_mask_width
			&& threadIdx.y >= blockDim.y -half_mask_width)
	{
		s_pix[threadIdx.y -(blockDim.y - half_mask_width)][threadIdx.x-(blockDim.x-half_mask_width)]
		      = (halo_index_x < 0 || halo_index_y < 0)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();

	//load the upper middle area
	halo_index_x = id_x;
	halo_index_y = (blockIdx.y-1)*blockDim.y + threadIdx.y;

	if(threadIdx.y >= blockDim.y -half_mask_width)
	{
		s_pix[threadIdx.y -(blockDim.y - half_mask_width)][half_mask_width + threadIdx.x]
		      = (halo_index_x >= width || halo_index_y < 0)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();

	//load the upper right area
	halo_index_x = (blockIdx.x+1)*blockDim.x + threadIdx.x;
	halo_index_y = (blockIdx.y-1)*blockDim.y + threadIdx.y;

	if(threadIdx.x < half_mask_width && threadIdx.y >= blockDim.y -half_mask_width)
	{
		s_pix[threadIdx.y -(blockDim.y - half_mask_width)][half_mask_width + blockDim.x + threadIdx.x]
		      = (halo_index_x >= width || halo_index_y < 0)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();

	//load the middle left area
	halo_index_x = (blockIdx.x-1)*blockDim.x + threadIdx.x;
	halo_index_y = id_y;

	if(threadIdx.x >= blockDim.x - half_mask_width)
	{
		s_pix[half_mask_width+threadIdx.y][threadIdx.x-(blockDim.x-half_mask_width)]
		      = (halo_index_x < 0 || halo_index_y >= height)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();



	//load the middle right
	halo_index_x = (blockIdx.x+1)*blockDim.x + threadIdx.x;
	halo_index_y = id_y;

	if(threadIdx.x < half_mask_width)
	{
		s_pix[half_mask_width + threadIdx.y][half_mask_width + blockDim.x + threadIdx.x]
		      = (halo_index_x >= width || halo_index_y >= height)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();

	//load the bottom left
	halo_index_x = (blockIdx.x-1)*blockDim.x + threadIdx.x;
	halo_index_y = (blockIdx.y+1)*blockDim.y + threadIdx.y;

	if(threadIdx.x >= blockDim.x - half_mask_width
			&& threadIdx.y < half_mask_width)
	{
		s_pix[half_mask_width + blockDim.y + threadIdx.y][threadIdx.x-(blockDim.x-half_mask_width)]
		      = (halo_index_x < 0 || halo_index_y >= height)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();

	//load the bottom middle
	halo_index_x = id_x;
	halo_index_y = (blockIdx.y+1)*blockDim.y + threadIdx.y;

	if(threadIdx.y < half_mask_width)
	{
		s_pix[half_mask_width + blockDim.y + threadIdx.y][half_mask_width + threadIdx.x]
		      = (halo_index_y >= height || halo_index_x >= width)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();

	//load the bottom right
	halo_index_x = (blockIdx.x+1)*blockDim.x + threadIdx.x;
	halo_index_y = (blockIdx.y+1)*blockDim.y + threadIdx.y;

	if(threadIdx.x < half_mask_width && threadIdx.y < half_mask_width)
	{
		s_pix[half_mask_width + blockDim.y + threadIdx.y][half_mask_width + blockDim.x + threadIdx.x]
		      = (halo_index_y >= height || halo_index_x >= width)? 0:d_pixArrayIn[halo_index_y*width + halo_index_x];
	}
	__syncthreads();

	//load the middle middle area
	halo_index_x = id_x;
	halo_index_y = id_y;

	s_pix[half_mask_width + threadIdx.y][half_mask_width + threadIdx.x]
	      = (halo_index_x < width && halo_index_y < height)?d_pixArrayIn[halo_index_y*width + halo_index_x]:0;
	__syncthreads();

	if(id_x < width && id_y < height)
	{
		float pvalue = 0;
		for(int i=0; i<mask_width; i++)
			for(int j=0; j<mask_width; j++)
				pvalue += s_pix[threadIdx.y+i][threadIdx.x+j] * c_mask[i*mask_width+j];

		d_pixArrayOut[id_y*width + id_x] = pvalue;
	}
}


extern "C" void cudaFilter(float *h_mask, unsigned int mask_width)
{
	float time;
	cudaEventRecord(start, 0);

	if(mask_width != MASK_WIDTH)
	{
		printf("convolve mask width configuration does not match the macro!\n");
		exit(-1);
	}

	if(mask_width/2 >= SQUARE_BLK_WIDTH)
	{
		printf("convolve mask width is too large, "
				"please use larger convolve kernel block size!\n");
		exit(-1);
	}

	checkCudaErrors(cudaMemcpyToSymbol(c_mask, h_mask, mask_width*mask_width*sizeof(float)));

	dim3 grid((width - 1) / SQUARE_BLK_WIDTH + 1,
			(height - 1) / SQUARE_BLK_HEIGHT + 1);
	dim3 block(SQUARE_BLK_WIDTH, SQUARE_BLK_HEIGHT);

	convolve_kernel<<<grid, block>>>(d_pixelArray, d_cdPixArray, width, height, mask_width);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

#ifdef DEBUG_CUDA
	printf("Time counsumed by convolve filter is: %f\n", time);
#endif

}


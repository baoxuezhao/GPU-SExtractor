#include <stdio.h>
#include <cuda.h>
#include <limits.h>
#include <iostream>
#include "../util5/cudpp.h"
#include "../util5/helper_cuda.h"

#include "cudatypes.h"

__global__ void compute_amp_alpha_kernel(
		float 	*d_a,
		float	*d_b,
		float	*d_fdflux,
		float	*d_abcor,
		float	*d_dthresh,
		unsigned int 	*d_fdnpix,
		double	*d_amp, 	//out
		double	*d_alpha, 	//out
		double	clean_param,
		unsigned int  numObj) {

	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid < numObj) {
		float unitareain = PI*d_a[tid]*d_b[tid];
		double ampin = d_amp[tid] = d_fdflux[tid]/(2*unitareain*d_abcor[tid]);
		d_alpha[tid] = (pow(ampin/d_dthresh[tid], 1.0/clean_param)-1)*unitareain/d_fdnpix[tid];
	}
}

__global__ void crosscompare_kernel(double *d_mx,
		double 	*d_my,
		float	*d_a,
		float	*d_fdflux,
		float	*d_cxx,
		float	*d_cyy,
		float	*d_cxy,
		float	*d_mthresh,
		double	*d_amp,
		double	*d_alpha,
		unsigned int *d_masterIndex,
		double	clean_param,
		unsigned int numObj) {

	//thread id within the whole grid
	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	//__shared__ double mx[2*SQUARE_BLK_WIDTH]; 	//256 bytes
	//__shared__ double my[2*SQUARE_BLK_WIDTH]; 	//256 bytes
	//__shared__ float a[2*SQUARE_BLK_WIDTH];		//128 bytes

	if(tid_x < numObj && tid_y < numObj) {

		//thread id within a block
		int id_x = threadIdx.x;
		int id_y = threadIdx.y;
		/*
		//load data into shared memory
		if(threadIdx.y == 0)
		{
			mx[threadIdx.x] = d_mx[tid_x];
			mx[threadIdx.x+blockDim.x] = d_mx[id_x+blockIdx.y*blockDim.y];

			my[threadIdx.x] = d_my[tid_x];
			my[threadIdx.x+blockDim.x] = d_my[id_x+blockIdx.y*blockDim.y];

			a[threadIdx.x] = d_a[tid_x];
			a[threadIdx.x+blockDim.x] = d_a[id_x+blockIdx.y*blockDim.y];
		}
		__syncthreads();*/

		if(tid_x == tid_y)
			return;

		//float dx = mx[threadIdx.x] - mx[blockDim.x+threadIdx.y];
		//float dy = my[threadIdx.x] - my[blockDim.x+threadIdx.y];
		//float rlim = a[threadIdx.x]+ a[blockDim.x+threadIdx.y];

		float dx = d_mx[tid_x] - d_mx[tid_y];
		float dy = d_my[tid_x] - d_my[tid_y];
		float rlim = d_a[tid_x] + d_a[tid_y];

		double val;

		rlim *= rlim;

		//row(obj), column(objin)
		if (dx*dx+dy*dy<rlim*CLEAN_ZONE*CLEAN_ZONE) {

			if(d_fdflux[tid_x] < d_fdflux[tid_y]) {

				val = 1+d_alpha[tid_y]*(d_cxx[tid_y]*dx*dx+d_cyy[tid_y]*dy*dy+d_cxy[tid_y]*dx*dy);

				//if (val>1.0 && ((float)(val<1e10?d_amp[tid_y]*pow(val,-clean_param) : 0.0) > d_mthresh[tid_x]))
				if (val>1.0 && ((float)(val<1e10?d_amp[tid_y]*(1.0/val):0.0) > d_mthresh[tid_x]))
				{
					//atomicCAS(d_masterIndex+tid_x, -1, tid_y);
					d_masterIndex[tid_x] = tid_y; //time consuming
					//atomicAdd(d_masterIndex+tid_x, 1); //time consuming
				}
			}
			/*
			else
			{
				val = 1+d_alpha[tid_x]*(d_cxx[tid_x]*dx*dx+d_cyy[tid_x]*dy*dy+d_cxy[tid_x]*dx*dy);
				if (val>1.0 && ((float)(val<1e10?d_amp[tid_x]*pow(val,-clean_param) : 0.0) > d_mthresh[tid_y]))
				{
					//d_masterIndex[tid_y] = tid_x;
					//atomicAdd(d_masterIndex+tid_y, 1);
				}
			}*/
		}
	}
}

/**
 * shape: x:numobj, y:3000
 */
__global__ void crosscompare_kernel1(
		double  *d_mx,
		double 	*d_my,
		float	*d_a,
		float	*d_fdflux,
		float	*d_cxx,
		float	*d_cyy,
		float	*d_cxy,
		float	*d_mthresh,
		double	*d_amp,
		double	*d_alpha,
		unsigned int *d_masterIndex,
		unsigned int numObj) {

	//thread id within the whole grid
	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	int neighbour = tid_x + tid_y + 1;

	if(tid_x < numObj && neighbour < numObj) {

		float dx = d_mx[tid_x] - d_mx[neighbour];
		float dy = d_my[tid_x] - d_my[neighbour];
		float rlim = d_a[tid_x] + d_a[neighbour];

		double val;

		rlim *= rlim;

		//row(obj), column(objin)
		if (dx*dx+dy*dy<rlim*CLEAN_ZONE*CLEAN_ZONE) {

			if(d_fdflux[tid_x] < d_fdflux[neighbour]) {

				val = 1+d_alpha[neighbour]*(d_cxx[neighbour]*dx*dx
						+d_cyy[neighbour]*dy*dy+d_cxy[neighbour]*dx*dy);
				if (val>1.0 && ((float)(val<1e10?d_amp[neighbour]*(1.0/val):0.0) > d_mthresh[tid_x]))
				{
					//d_masterIndex[tid_x] = neighbour+10; //time consuming
					atomicMin(d_masterIndex+tid_x, neighbour);
				}
			}
			else
			{
				val = 1+d_alpha[tid_x]*(d_cxx[tid_x]*dx*dx+d_cyy[tid_x]*dy*dy+d_cxy[tid_x]*dx*dy);
				if (val>1.0 && ((float)(val<1e10?d_amp[tid_x]*(1.0/val) : 0.0) > d_mthresh[neighbour]))
				{
					//d_masterIndex[neighbour] = tid_x+10;
					atomicMin(d_masterIndex+neighbour, tid_x);
				}
			}
		}
	}
}


__global__ void mergeobject_kernel(
		unsigned int *d_masterIndex,
		unsigned int *d_sortedObjIndex,
		unsigned int *d_fdnpix,
		unsigned int *d_dnpix,
		float	*d_fdflux,
		float	*d_dflux,
		float	*d_flux,
		float	*d_fluxerr,
		float	*d_fdpeak,
		float	*d_dpeak,
		float	*d_peak,
		unsigned int *d_peakx,
		unsigned int *d_peaky,
		unsigned int *d_xmin,
		unsigned int *d_xmax,
		unsigned int *d_ymin,
		unsigned int *d_ymax,
		unsigned int *d_flag,
		unsigned int numObj) {

	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= numObj || d_masterIndex[tid] >= numObj)
		return;

	if(tid >0 && d_masterIndex[tid] == d_masterIndex[tid-1])
		return;

	int idx = tid;
	int masterIdx =d_masterIndex[idx];
	int slaveIdx;
	while(d_masterIndex[idx] == masterIdx)
	{
		slaveIdx = d_sortedObjIndex[idx];

		d_fdnpix[masterIdx] += d_fdnpix[slaveIdx];
		d_dnpix[masterIdx]  += d_dnpix[slaveIdx];
		d_fdflux[masterIdx] += d_fdflux[slaveIdx];
		d_dflux[masterIdx]  += d_dflux[slaveIdx];
		d_flux[masterIdx]   += d_flux[slaveIdx];
		d_fluxerr[masterIdx]+= d_fluxerr[slaveIdx];

		if (d_fdpeak[slaveIdx]>d_fdpeak[masterIdx])
		{
			d_fdpeak[masterIdx] = d_fdpeak[slaveIdx];
			d_peakx[masterIdx]  = d_peakx[slaveIdx];
			d_peaky[masterIdx]  = d_peaky[slaveIdx];
		}
		if (d_dpeak[slaveIdx] > d_dpeak[masterIdx])
			d_dpeak[masterIdx] = d_dpeak[slaveIdx];

		if (d_peak[slaveIdx] > d_peak[masterIdx])
			d_peak[masterIdx] = d_peak[slaveIdx];

		if (d_xmin[slaveIdx] < d_xmin[masterIdx])
			d_xmin[masterIdx] = d_xmin[slaveIdx];
		if (d_xmax[slaveIdx] > d_xmax[masterIdx])
			d_xmax[masterIdx] = d_xmax[slaveIdx];

		if (d_ymin[slaveIdx] < d_ymin[masterIdx])
			d_ymin[masterIdx] = d_ymin[slaveIdx];
		if (d_ymax[slaveIdx] > d_ymax[masterIdx])
			d_ymax[masterIdx] = d_ymax[slaveIdx];

		d_flag[masterIdx] |= (d_flag[slaveIdx] & (~(OBJ_MERGED|OBJ_CROWDED)));

		//prefs.nimaisoflag = 0 currently
		//mergeflags(objmaster, objslave);

		idx++;
	}

	//fdnpix, dnpix
	//fdflux, dflux, flux, fluxerr

	//fdpeak, peakx, peaky
	//dpeak, peak
	//xmin, xmax, ymin, ymax
	//flag
}


extern "C" void run_clean(int *masterIndex,
		double clean_param,
		int clean_stacksize,
		int numobj)
{

	//outputfortest(numobj);

	dim3 grid((numobj-1)/SQUARE_BLK_HEIGHT + 1, (clean_stacksize-1)/SQUARE_BLK_WIDTH + 1);
	dim3 block(SQUARE_BLK_WIDTH, SQUARE_BLK_HEIGHT);

	int grid1 = (numobj-1)/SQUARE_BLK_WIDTH + 1;
	int block1= SQUARE_BLK_WIDTH;

	double *d_amp, 	*d_alpha;
	unsigned int 	*d_sortedmasterIdx;

	float time;
	cudaEventRecord(start, 0);

	cudaMalloc((void**)&d_amp,		numobj*sizeof(double));
	cudaMalloc((void**)&d_alpha,	numobj*sizeof(double));
	cudaMalloc((void**)&d_masterIndex, numobj*sizeof(unsigned int));
	cudaMalloc((void**)&d_sortedmasterIdx, numobj*sizeof(unsigned int));
	cudaMemset(d_masterIndex, -1, 	numobj*sizeof(unsigned int));

	compute_amp_alpha_kernel<<<grid1, block1>>>(
			d_a,
			d_b,
			d_fdflux,
			d_abcor,
			d_dthresh,
			d_fdnpix,
			d_amp, 		//out
			d_alpha, 	//out
			clean_param,
			numobj);

	crosscompare_kernel1<<<grid, block>>>(d_mx,
			d_my,
			d_a,
			d_fdflux,
			d_cxx,
			d_cyy,
			d_cxy,
			d_mthresh,
			d_amp,
			d_alpha,
			d_masterIndex,
			numobj);

	/*
	cudaMemcpy(d_sortedmasterIdx, d_masterIndex, numobj*sizeof(unsigned int), cudaMemcpyDeviceToDevice);

	//d_cuttedObjLabelArray is initialized in analysis_kernel
	CUDPPResult res;

	//sort the cutted object by root label
	//CUDPPConfiguration sort_config;
	config.datatype = CUDPP_UINT;
	config.algorithm = CUDPP_SORT_RADIX;
	config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;

	//create the cudpp sort plan
	//CUDPPHandle sortplan = 0;
	res = cudppPlan(theCudpp, &sortplan, config, numobj, 1, 0);
	if (CUDPP_SUCCESS != res)
	{
		printf("Error creating CUDPP sort Plan in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	// Run the sort
	res = cudppSort(sortplan, d_sortedmasterIdx, d_cuttedObjLabelArray, numobj);
	if (CUDPP_SUCCESS != res)
	{
		printf("Error in cudppSort() for sorting d_masterIndex in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	// Destroy the sort plan
	res = cudppDestroyPlan(sortplan);
	if (CUDPP_SUCCESS != res)
	{
		printf("Error destroying CUDPP sort Plan in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	 mergeobject_kernel<<<grid1, block1>>>(
			 	d_sortedmasterIdx,
				d_cuttedObjLabelArray,
				d_fdnpix,
				d_dnpix,
				d_fdflux,
				d_dflux,
				d_flux,
				d_fluxerr,
				d_fdpeak,
				d_dpeak,
				d_peak,
				d_peakx,
				d_peaky,
				d_xmin,
				d_xmax,
				d_ymin,
				d_ymax,
				d_flag,
				numobj);
	*/

	//int count = 0;
	cudaMemcpy(masterIndex, d_masterIndex, numobj*sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

	printf("Time counsumed by cuda clean is: %f\n", time);

#ifdef DEBUG_CUDA
	printf("Time counsumed by cuda clean is: %f\n", time);
#endif

	cudaFree(d_amp);
	cudaFree(d_alpha);
	cudaFree(d_sortedmasterIdx);
}


#include <stdio.h>
#include <cuda.h>
#include <limits.h>
#include <math.h>
#include "../util5/cudpp.h"
#include "../util5/helper_cuda.h"

#include "cudatypes.h"

#define INBLOCK_IDX(i,j) 	(threadIdx.x+(i)) + (threadIdx.y+(j))*blockDim.x
#define BLOCK_IDX			blockIdx.x + blockIdx.y*gridDim.x


#define CUDA_CHECK_RETURN(value) {											\
		cudaError_t _m_cudaStat = value;										\
		if (_m_cudaStat != cudaSuccess) {										\
			fprintf(stderr, "Error %s at line %d in file %s\n",					\
					cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
					exit(1);															\
		} }

void get_mean_sigma(float *d_mean1,
		float *d_sigma1,
		float *backmean,
		float *backsig,
		int nb);
/*
Median using Heapsort algorithm (for float arrays) (based on Num.Rec algo.).
Warning: changes the order of data!
 */
__device__ float hmedian(float *ra, int n)
{
	int		l, j, ir, i;
	float	rra;


	if (n<2)
		return *ra;
	ra--;
	for (l = ((ir=n)>>1)+1;;)
	{
		if (l>1)
			rra = ra[--l];
		else
		{
			rra = ra[ir];
			ra[ir] = ra[1];
			if (--ir == 1)
			{
				ra[1] = rra;
				return n&1? ra[n/2+1] : (float)((ra[n/2]+ra[n/2+1])/2.0);
			}
		}
		for (j = (i=l)<<1; j <= ir;)
		{
			if (j < ir && ra[j] < ra[j+1])
				++j;
			if (rra < ra[j])
			{
				ra[i] = ra[j];
				j += (i=j);
			}
			else
				j = ir + 1;
		}
		ra[i] = rra;
	}

	/* (the 'return' is inside the loop!!) */
}

/**
 * thread size: nmesh_x * nmesh_y
 * we assume that the mesh_size >= blockDim.x(y)
 */
__global__ void backstat_kernel(float *d_imgbuf,
		int 	*d_npix,
		int 	*d_nlevels,
		float 	*d_mean,
		float 	*d_sigma,
		float 	*d_lcut,
		float 	*d_hcut,
		float	*d_qscale,
		float	*d_qzero,
		int		*d_histo,
		int		width,
		int		height,
		int		back_size)
{

	int i,j;
	int nx = (back_size-1)/blockDim.x + 1;
	int ny = (back_size-1)/blockDim.y + 1;

	int innermesh_x;
	int innermesh_y;

	float  pixf, qscale, cste;
	double mean1, sigma1, pix;
	double step = sqrt(2/PI)*QUANTIF_NSIGMA/QUANTIF_AMIN;

	//size: (2*blockDim.x*blockDim.y)*sizeof(float) + 3*sizeof(int)
	extern __shared__ char shared_space[];

	//store pixels(back_size x back_size)
	float *pix_t 	= (float*)shared_space;

	//store mean values (blockDim.x * blockDim.y)
	double *mean_t  = (double*)&shared_space[back_size*back_size*sizeof(float)];

	//store sigma values(blockDim.x * blockDim.y)
	double *sigma_t = (double*)&shared_space[(back_size*back_size)*sizeof(float)
	                                         +blockDim.x*blockDim.y*sizeof(double)];

	//store count(good pixel) (blockDim.x * blockDim.y)
	int   *npix_t = (int*)&shared_space[(back_size*back_size)*sizeof(float)
	                                    +2*blockDim.x*blockDim.y*sizeof(double)];

	//store count(all pixel) (blockDim.x * blockDim.y)
	int   *npix_tt = (int*)&shared_space[(back_size*back_size)*sizeof(float)
	                                    +2*blockDim.x*blockDim.y*sizeof(double)
										+blockDim.x*blockDim.y*sizeof(int)];
	//1
	float *lcut_t = (float*)&shared_space[(back_size*back_size)*sizeof(float)
	                                    +2*blockDim.x*blockDim.y*sizeof(double)
										+2*blockDim.x*blockDim.y*sizeof(int)];
	//1
	float *hcut_t = (float*)&shared_space[(back_size*back_size)*sizeof(float)
	                                    +2*blockDim.x*blockDim.y*sizeof(double)
										+2*blockDim.x*blockDim.y*sizeof(int)+sizeof(float)];
	//init to 0
	mean_t [INBLOCK_IDX(0,0)] = 0.0;
	sigma_t[INBLOCK_IDX(0,0)] = 0.0;
	npix_t [INBLOCK_IDX(0,0)] = 0;
	npix_tt[INBLOCK_IDX(0,0)] = 0;

	for(int i=0; i<nx; ++i) {
		for(int j=0; j<ny; ++j) {

			innermesh_x = (blockIdx.x*back_size) + (i*blockDim.x) + threadIdx.x;
			innermesh_y = (blockIdx.y*back_size) + (j*blockDim.y) + threadIdx.y;

			if((i*blockDim.x) + threadIdx.x < back_size
				&& (j*blockDim.y) + threadIdx.y < back_size
				&& innermesh_x < width && innermesh_y < height)
			{
				pix = pix_t[(i*blockDim.x) + threadIdx.x + ((j*blockDim.y) + threadIdx.y)*back_size] =
						d_imgbuf[innermesh_x + innermesh_y*width];
				if(pix > -BIG) {
					mean_t [INBLOCK_IDX(0,0)] += pix;
					sigma_t[INBLOCK_IDX(0,0)] += (pix*pix);
					npix_t [INBLOCK_IDX(0,0)]++;
				}
				npix_tt[INBLOCK_IDX(0,0)]++;
			}
		}
	}
	__syncthreads();

	for(i=blockDim.x/2; i>0; i/=2)
	{
		if(threadIdx.x < i)
		{
			mean_t [INBLOCK_IDX(0,0)] += mean_t [INBLOCK_IDX(i,0)];
			sigma_t[INBLOCK_IDX(0,0)] += sigma_t[INBLOCK_IDX(i,0)];
			npix_t [INBLOCK_IDX(0,0)] += npix_t [INBLOCK_IDX(i,0)];
			npix_tt [INBLOCK_IDX(0,0)] += npix_tt [INBLOCK_IDX(i,0)];
		}
		__syncthreads();
	}

	if(threadIdx.x == 0) {
		for(j=blockDim.y/2; j>0; j/=2)
		{
			if(threadIdx.y < j) {
				mean_t [INBLOCK_IDX(0,0)] += mean_t [INBLOCK_IDX(0,j)];
				sigma_t[INBLOCK_IDX(0,0)] += sigma_t[INBLOCK_IDX(0,j)];
				npix_t [INBLOCK_IDX(0,0)] += npix_t [INBLOCK_IDX(0,j)];
				npix_tt[INBLOCK_IDX(0,0)] += npix_tt[INBLOCK_IDX(0,j)];
			}
			__syncthreads();
		}
	}

	if(threadIdx.x == 0 && threadIdx.y == 0)
	{
		if ((float)npix_t[0] < (float)(npix_tt[0]*BACK_MINGOODFRAC))
		{
			d_mean[BLOCK_IDX] = d_sigma[BLOCK_IDX] =  -BIG;
			//????????????
			return;
		}
		mean1 = mean_t[0] /(double)npix_t[0];
		sigma1= (sigma_t[1]=sigma_t[0]/npix_t[0]-mean1*mean1)>0.0?sqrt(sigma_t[1]):0.0;

		lcut_t[0] = (mean1 - 2.0*sigma1);
		hcut_t[0] = (mean1 + 2.0*sigma1);
		__syncthreads();
	}

	npix_t[INBLOCK_IDX(0,0)] = 0;
	mean_t[INBLOCK_IDX(0,0)] = 0.0;
	sigma_t[INBLOCK_IDX(0,0)] = 0.0;

	__syncthreads();

	for(int i=0; i<nx; ++i) {
		for(int j=0; j<ny; ++j) {

			innermesh_x = (blockIdx.x*back_size) + (i*blockDim.x) + threadIdx.x;
			innermesh_y = (blockIdx.y*back_size) + (j*blockDim.y) + threadIdx.y;

			if((i*blockDim.x) + threadIdx.x < back_size
				&& (j*blockDim.y) + threadIdx.y < back_size
				&& innermesh_x < width && innermesh_y < height)
			{
				pix = pix_t[(i*blockDim.x) + threadIdx.x + ((j*blockDim.y) + threadIdx.y)*back_size];
				if(pix<=hcut_t[0] && pix>=lcut_t[0]) {
					mean_t [INBLOCK_IDX(0,0)] += pix;
					sigma_t[INBLOCK_IDX(0,0)] += (pix*pix);
					npix_t [INBLOCK_IDX(0,0)]++;
				}
			}
		}
	}

	for(i=blockDim.x/2; i>0; i/=2)
	{
		if(threadIdx.x < i)
		{
			mean_t[INBLOCK_IDX(0,0)]  += mean_t[INBLOCK_IDX(i,0)];
			sigma_t[INBLOCK_IDX(0,0)] += sigma_t[INBLOCK_IDX(i,0)];
			npix_t[INBLOCK_IDX(0,0)]  += npix_t[INBLOCK_IDX(i,0)];
		}
	}
	__syncthreads();

	if(threadIdx.x == 0) {
		for(j=blockDim.y/2; j>0; j/=2)
		{
			if(threadIdx.y < j)
			{
				mean_t[INBLOCK_IDX(0,0)]  += mean_t[INBLOCK_IDX(0,j)];
				sigma_t[INBLOCK_IDX(0,0)] += sigma_t[INBLOCK_IDX(0,j)];
				npix_t[INBLOCK_IDX(0,0)]  += npix_t[INBLOCK_IDX(0,j)];
			}
		}
		__syncthreads();
	}

	if(threadIdx.x == 0 && threadIdx.y == 0)
	{
		d_npix[BLOCK_IDX] = npix_t[0];
		d_mean[BLOCK_IDX] = mean1 = (mean_t[0]/(double)npix_t[0]);
		d_sigma[BLOCK_IDX] = sigma1 =
				(sigma_t[1]=sigma_t[0]/(double)npix_t[0]-mean1*mean1)>0.0?sqrt(sigma_t[1]):0.0;

		if((d_nlevels[BLOCK_IDX] = (int)(step*npix_t[0]+1))>QUANTIF_NMAXLEVELS)
			d_nlevels[BLOCK_IDX] = QUANTIF_NMAXLEVELS;

		d_qscale[BLOCK_IDX] = qscale =
				sigma1>0.0? 2*QUANTIF_NSIGMA*sigma1/d_nlevels[BLOCK_IDX] : 1.0;
		d_qzero[BLOCK_IDX] = mean1 - QUANTIF_NSIGMA*sigma1;

		d_lcut[BLOCK_IDX] = lcut_t[0];
		d_hcut[BLOCK_IDX] = hcut_t[0];
	}
	__syncthreads();

	//backhisto
	//int[0] to int[4095] takes up all the shared memory(48k) except pix_t
	int *histo = (int*)&shared_space[/*2*/back_size*back_size*sizeof(float)];

	for(int k=0; k<QUANTIF_NMAXLEVELS/(blockDim.x*blockDim.y); k++)
	{
		if(k*(blockDim.x*blockDim.y) + INBLOCK_IDX(0,0) < QUANTIF_NMAXLEVELS)
		histo[k*(blockDim.x*blockDim.y) + INBLOCK_IDX(0,0)] = 0;
	}
	__syncthreads();

	//////////////////////////
	//int bindex = blockIdx.x + (blockIdx.y*gridDim.x);  //index of blocks
	//int offset = bindex*QUANTIF_NMAXLEVELS; //offset to write in the d_histo

	if(d_mean[BLOCK_IDX] > -BIG) {

		//memset(histo, 0, QUANTIF_NMAXLEVELS * sizeof(int));

		int bin = 0;
		qscale = d_qscale[BLOCK_IDX];
		cste = 0.499999 - d_qzero[BLOCK_IDX]/qscale;

		for(int i=0; i<nx; ++i) {
			for(int j=0; j<ny; ++j) {

				innermesh_x = (blockIdx.x*back_size) + (i*blockDim.x) + threadIdx.x;
				innermesh_y = (blockIdx.y*back_size) + (j*blockDim.y) + threadIdx.y;

				if((i*blockDim.x) + threadIdx.x < back_size
					&& (j*blockDim.y) + threadIdx.y < back_size
					&& innermesh_x < width && innermesh_y < height)
				{
					pixf = pix_t[(i*blockDim.x) + threadIdx.x + ((j*blockDim.y) + threadIdx.y)*back_size];
					bin = (int)(pixf/qscale + cste);

					if(bin >= 0 && bin < d_nlevels[BLOCK_IDX])
					{
						atomicAdd(histo+bin, 1);
					}
				}
				__syncthreads();
			}
		}
	}

	__syncthreads();

	//copy histo result back to device memory
	int bindex = BLOCK_IDX;  //index of blocks
	int offset = bindex*QUANTIF_NMAXLEVELS; //offset to write in the d_histo
	int nvalue= QUANTIF_NMAXLEVELS / (blockDim.x*blockDim.y); //num of values each thread need to write

	for(int k=0; k<nvalue; k++)
	{
		if(k*(blockDim.x*blockDim.y) + INBLOCK_IDX(0,0) < QUANTIF_NMAXLEVELS)
		d_histo[offset + k*(blockDim.x*blockDim.y) + INBLOCK_IDX(0,0)]
		        = histo[k*(blockDim.x*blockDim.y) + INBLOCK_IDX(0,0)];
	}


	//////////////////////////////////////////////////

}

__global__ void backguess_kernel(int *d_nlevels,
		float 	*d_mean,
		float 	*d_sigma,
		float 	*d_qscale,
		float 	*d_qzero,
		int		*d_histo,
		//float	*d_outmean,
		//float	*d_outsigma,
		int		nx_mesh,
		int		ny_mesh) {

	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	if(tid_x >= nx_mesh || tid_y >= ny_mesh)
		return;

	const int tid = tid_y * nx_mesh + tid_x;

	/* Leave here if the mesh is already classified as `bad' */
	if (d_mean[tid] <= -BIG)
	{
		d_mean[tid] = d_sigma[tid] = -BIG;
		return;
	}

	int		*histo, *hilow, *hihigh, *histot;
	unsigned long lowsum, highsum, sum;
	double	ftemp, mea, sig, sig1, med, dpix;
	int		i, n, lcut,hcut, nlevelsm1, pix;

	histo = d_histo + tid*QUANTIF_NMAXLEVELS;
	hcut = nlevelsm1 = d_nlevels[tid]-1;
	lcut = 0;

	sig = 10.0*nlevelsm1;
	sig1 = 1.0;
	mea = med = d_mean[tid];

	for (n=100; n-- && (sig>=0.1) && (fabs(sig/sig1-1.0)>EPS);)
	{
		sig1 = sig;
		sum = mea = sig = 0.0;
		lowsum = highsum = 0;
		histot = hilow = histo+lcut;
		hihigh = histo+hcut;

		for (i=lcut; i<=hcut; i++)
		//for (i=hcut; i>=lcut; i--)
		{
			if (lowsum<highsum)
				lowsum += *(hilow++);
			else
				highsum += *(hihigh--);

			sum += (pix = *(histot++));
			mea += (dpix = (double)pix*i);
			sig += dpix*i;
		}

		med = hihigh>=histo?
				((hihigh-histo)+0.5+((double)highsum-lowsum)/(2.0*(*hilow>*hihigh?
						*hilow:*hihigh))) : 0.0;

		//if(tid == 20)
		//	printf("%d, %f, %f, %f, %d \n", sum, mea, sig, med, hcut-lcut);

		if (sum)
		{
			mea /= (double)sum;
			sig = sig/sum - mea*mea;
		}
		sig = sig>0.0?sqrt(sig):0.0;
		lcut = (ftemp=med-3.0*sig)>0.0 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5):0;
		hcut = (ftemp=med+3.0*sig)<nlevelsm1 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5)
				: nlevelsm1;
	}

	if(abs(sig)>0.0)
	{
		if(fabs(d_sigma[tid]/(sig*d_qscale[tid])-1) < 0.0) {

			d_mean[tid] = d_qzero[tid]+mea*d_qscale[tid];
		}
		else
		{
			if(fabs((mea-med)/sig)< 0.3)
				d_mean[tid] = d_qzero[tid]+(2.5*med-1.5*mea)*d_qscale[tid];
			else
				d_mean[tid] = d_qzero[tid]+med*d_qscale[tid];
		}
	}else
	{
		d_mean[tid] = d_qzero[tid]+mea*d_qscale[tid];
	}

	d_sigma[tid] = sig*d_qscale[tid];
}


//not used currently
__global__ void backguess_kernel1(
		int 	*d_nlevels,
		float 	*d_mean,
		float 	*d_sigma,
		float 	*d_qscale,
		float 	*d_qzero,
		int		*d_histo,
		int		*d_lowsum,
		int		*d_highsum,
		int		num_mesh,
		int		histo_size)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= histo_size)
		return;

	//int histoIdx = tid/QUANTIF_NMAXLEVELS;
	int nx = (QUANTIF_NMAXLEVELS-1)/blockDim.x+1;

	//if(tid_x >= nx_mesh || tid_y >= ny_mesh)
	//	return;

	//const int bid = blockIdx.y * nx_mesh + blockIdx.x;

	if(threadIdx.x == 0 && threadIdx.y == 0)
	{
		/* Leave here if the mesh is already classified as `bad' */
		if (d_mean[BLOCK_IDX] <= -BIG)
		{
			d_mean[BLOCK_IDX] = d_sigma[BLOCK_IDX] = -BIG;
			return;
		}
	}

	__syncthreads();

	if(d_mean[BLOCK_IDX] == -BIG && d_sigma[BLOCK_IDX] == -BIG)
		return;

	int *histo 	= d_histo + BLOCK_IDX*QUANTIF_NMAXLEVELS;
	int *lowsum = d_lowsum + BLOCK_IDX*QUANTIF_NMAXLEVELS;
	int *highsum= d_highsum + BLOCK_IDX*QUANTIF_NMAXLEVELS;

	extern __shared__ char shared_space[];

	//16k
	int	*histo_t = (int*)shared_space;

	//4k
	int *hilowdiff_t = (int*)&shared_space[QUANTIF_NMAXLEVELS*sizeof(int)];

	//4k
	int *diffpos_t   = (int*)&shared_space[QUANTIF_NMAXLEVELS*sizeof(int)
									+ blockDim.x*blockDim.y*sizeof(int)];
	//4k
	int *sum_t = (int*)&shared_space[QUANTIF_NMAXLEVELS*sizeof(int)
									+ 2*blockDim.x*blockDim.y*sizeof(int)];
	//8k
	double *mea_t = (double*)&shared_space[QUANTIF_NMAXLEVELS*sizeof(int)
									+ 3*blockDim.x*blockDim.y*sizeof(int)];
	//8k
	double *sig_t = (double*)&shared_space[QUANTIF_NMAXLEVELS*sizeof(int)
									+ 3*blockDim.x*blockDim.y*sizeof(int)
									+ blockDim.x*blockDim.y*sizeof(double)];
	//4 byte
	int	*l_cut_t = (int*)&shared_space[QUANTIF_NMAXLEVELS*sizeof(int)
									+ 3*blockDim.x*blockDim.y*sizeof(int)
									+ 2*blockDim.x*blockDim.y*sizeof(double)];
	//4 byte
	int *h_cut_t = (int*)&shared_space[QUANTIF_NMAXLEVELS*sizeof(int)
									+ 3*blockDim.x*blockDim.y*sizeof(int)
									+ 2*blockDim.x*blockDim.y*sizeof(double)
									+ sizeof(int)];

	int index, n, pix, diff, low, high;
	double sig, sig1, dpix;
	sig = 10.0*(d_nlevels[BLOCK_IDX]-1);
	sig1 = 1.0;

	for (n=100; n-- && (sig>=0.1) && (fabs(sig/sig1-1.0)>EPS);)
	{
		//init to 0
		hilowdiff_t [threadIdx.x] = INT_MAX;
		sum_t [threadIdx.x] = 0;
		mea_t [threadIdx.x] = 0.0;
		sig_t [threadIdx.x] = 0.0;

		if(threadIdx.x == 0 && threadIdx.y == 0)
		{
			l_cut_t[0] = 0;
			h_cut_t[0] = d_nlevels[BLOCK_IDX]-1;
		}
		__syncthreads();

		for(int i=0; i<nx; ++i)
		{
			index = (i*blockDim.x) + threadIdx.x;

			if(index >= l_cut_t[0] && index <= h_cut_t[0])
			{
				low = lowsum[index] - l_cut_t[0]?0:lowsum[l_cut_t[0]-1];
				high = highsum[index] -
						(h_cut_t[0] == d_nlevels[BLOCK_IDX]-1)?0:highsum[h_cut_t[0]+1];

				diff = abs(low-high);
				if(diff < hilowdiff_t[threadIdx.x])
				{
					diffpos_t[threadIdx.x] = index;
					hilowdiff_t[threadIdx.x] = diff;
				}

				sum_t[threadIdx.x] += (pix = histo[index]);
				mea_t[threadIdx.x] += (dpix = (double)pix*index);
				sig_t[threadIdx.x] += dpix*index;
			}
		}
		__syncthreads();

		//parallel reduction
		for(int i=blockDim.x/2; i>0; i/=2)
		{
			if(threadIdx.x < i)
			{
				if(hilowdiff_t[threadIdx.x] > hilowdiff_t[threadIdx.x+i])
				{
					hilowdiff_t[threadIdx.x] = hilowdiff_t[threadIdx.x+i];
					diffpos_t[threadIdx.x] = diffpos_t[threadIdx.x+i];
				}

				sum_t[threadIdx.x] += sum_t[threadIdx.x+i];
				mea_t[threadIdx.x] += mea_t[threadIdx.x+i];
				sig_t[threadIdx.x] += sig_t[threadIdx.x+i];
			}
			__syncthreads();
		}

	}

}

//not used currently
__global__ void initHistoSegment_kernel(
		int		*d_histosegment,
		int		histo_size)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= histo_size)
		return;

	d_histosegment[tid] = tid/QUANTIF_NMAXLEVELS;
}

//not used currently
__global__ void initMask_kernel(
		unsigned int *d_scanmask,
		int		*d_histosegment,
		int		histo_size)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= histo_size)
		return;

	if(tid == 0)
		d_scanmask[tid] = 1;
	else
	{
		if(d_histosegment[tid] > d_histosegment[tid-1])
			d_scanmask[tid] = 1;
		else
			d_scanmask[tid] = 0;
	}
}

//not used currently
__global__ void initdpix_kernel(
		int		*d_histo,
		double	*d_dpix,
		double	*d_dpixmi,
		int		*d_histosegment,
		int		*d_lcut,
		int		*d_startpos,
		int		histo_size)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= histo_size)
		return;

	int histoIdx = d_histosegment[tid];
	int offset = d_lcut[histoIdx] + (tid-d_startpos[histoIdx]);

	d_dpix[tid] = (double)d_histo[tid]*offset;
	d_dpixmi[tid] = d_dpix[tid]*offset;
}

//not used currently
__global__ void computeHighSum_kernel(
		int		*d_lowsum,
		int		*d_highsum,
		int		*d_hilowdiff,
		int		*d_lcut,
		int		*d_hcut,
		int		*d_histosegment,
		int		*d_startpos,
		int		histo_size)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= histo_size)
		return;

	int histoIdx = d_histosegment[tid];

	int endpos = d_startpos[histoIdx] + (d_hcut[histoIdx]-d_lcut[histoIdx]);

	int total = d_lowsum[endpos];

	d_highsum[tid] = total - d_lowsum[tid];

	d_hilowdiff[tid] = abs(d_highsum[tid] - d_lowsum[tid]);
}

//not used currently
__global__ void initBackguess_kernel(
		int		*d_lcut,
		int		*d_hcut,
		int		*d_startpos,
		int 	*d_nlevels,
		double 	*d_sig,
		double 	*d_sig1,
		int		num_mesh)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= num_mesh)
		return;

	int nlevelsm1;

	d_lcut[tid] = 0;
	d_hcut[tid] = nlevelsm1 = d_nlevels[tid]-1;;

	d_startpos[tid] = tid*QUANTIF_NMAXLEVELS;

	d_sig[tid] = 10.0 * nlevelsm1;
	d_sig1[tid] = 1.0;
}

//not used currently
void backguess_gpu(
		int 	*d_nlevels,
		float 	*d_mean,
		float	*d_sigma,
		float	*d_qscale,
		float	*d_qzero,
		int		*d_histo,
		int		nx_mesh,
		int		ny_mesh)
{
	int histo_size = nx_mesh * ny_mesh * QUANTIF_NMAXLEVELS;

	unsigned int 	*d_scanmask;
	int		*d_lowsum, *d_highsum, *d_hilowdiff;
	int		*d_lcut, *d_hcut, *d_startpos;
	int		*d_histosegment;
	double	*d_sig, *d_sig1;
	double 	*d_dpix, *d_dpixmi;

	double	*d_meansegment, *d_sigsegment;

	cudaMalloc((void**)&d_scanmask, histo_size*sizeof(unsigned int));
	cudaMalloc((void**)&d_lowsum, 	histo_size*sizeof(int));
	cudaMalloc((void**)&d_highsum, 	histo_size*sizeof(int));
	cudaMalloc((void**)&d_hilowdiff,histo_size*sizeof(int));
	cudaMalloc((void**)&d_histosegment, histo_size*sizeof(int));

	cudaMalloc((void**)&d_dpix, 	histo_size*sizeof(double));
	cudaMalloc((void**)&d_dpixmi, 	histo_size*sizeof(double));
	cudaMalloc((void**)&d_meansegment, 	histo_size*sizeof(double));
	cudaMalloc((void**)&d_sigsegment, 	histo_size*sizeof(double));

	cudaMemset(d_scanmask, 0, histo_size*sizeof(unsigned int));

	cudaMalloc((void**)&d_startpos, nx_mesh*ny_mesh*sizeof(int));
	cudaMalloc((void**)&d_lcut, nx_mesh*ny_mesh*sizeof(int));
	cudaMalloc((void**)&d_hcut, nx_mesh*ny_mesh*sizeof(int));
	cudaMalloc((void**)&d_sig,  nx_mesh*ny_mesh*sizeof(double));
	cudaMalloc((void**)&d_sig1, nx_mesh*ny_mesh*sizeof(double));

	int grid = (nx_mesh*ny_mesh-1)/MAX_THREADS_PER_BLK + 1;
	int block = MAX_THREADS_PER_BLK;

	initBackguess_kernel<<<grid, block>>>(
			d_lcut,
			d_hcut,
			d_startpos,
			d_nlevels,
			d_sig,
			d_sig1,
			nx_mesh*ny_mesh);

	grid = (histo_size-1)/MAX_THREADS_PER_BLK + 1;
	block = MAX_THREADS_PER_BLK;

	initHistoSegment_kernel<<<grid, block>>>(
			d_histosegment,
			histo_size);

	//for (n=100; n-- && (sig>=0.1) && (fabs(sig/sig1-1.0)>EPS);)

	int n = 1, go = 1;
	while(n && go)
	{
		//compute the sum, mea, sig, lowsum, highsum, lcut and hcut

		int grid1 = (histo_size-1)/MAX_THREADS_PER_BLK + 1;
		int block1 = MAX_THREADS_PER_BLK;

		initMask_kernel<<<grid1, block1>>>(
				d_scanmask,
				d_histosegment,
				histo_size);

		initdpix_kernel<<<grid1, block1>>>(
				d_histo,
				d_dpix,
				d_dpixmi,
				d_histosegment,
				d_lcut,
				d_startpos,
				histo_size);

		config.op = CUDPP_ADD;
		config.datatype = CUDPP_INT;
		config.algorithm = CUDPP_SEGMENTED_SCAN;
		config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE;

		CUDPPResult res = cudppPlan(theCudpp, &scanplan, config, histo_size, 1, 0);
		if (CUDPP_SUCCESS != res)
		{
			printf("Error in cudppSort() for sorting compacted labels in %s at line %d\n", __FILE__, __LINE__);
			exit(-1);
		}

		res = cudppSegmentedScan(scanplan,
				d_lowsum,
				d_histo,
				d_scanmask,
				histo_size);
		if (CUDPP_SUCCESS != res)
		{
			printf("Error in cudppSort() for sorting compacted labels in %s at line %d\n", __FILE__, __LINE__);
			exit(-1);
		}

		// Destroy the segment prefix scan plan
		res = cudppDestroyPlan(scanplan);
		if (CUDPP_SUCCESS != res)
		{
			printf("Error in cudppSort() for sorting compacted labels in %s at line %d\n", __FILE__, __LINE__);
			exit(-1);
		}

		computeHighSum_kernel<<<grid1, block1>>>(
				d_lowsum,
				d_highsum,
				d_hilowdiff,
				d_lcut,
				d_hcut,
				d_histosegment,
				d_startpos,
				histo_size);

		///////////////////////////////////
		config.op = CUDPP_ADD;
		config.datatype = CUDPP_DOUBLE;
		config.algorithm = CUDPP_SEGMENTED_SCAN;
		config.options = CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;

		res = cudppPlan(theCudpp, &scanplan, config, histo_size, 1, 0);
		if (CUDPP_SUCCESS != res)
		{
			printf("Error in cudppSort() for sorting compacted labels in %s at line %d\n", __FILE__, __LINE__);
			exit(-1);
		}

		//compute mean
		res = cudppSegmentedScan(scanplan,
				d_meansegment,
				d_dpix,
				d_scanmask,
				histo_size);
		if (CUDPP_SUCCESS != res)
		{
			printf("Error in cudppSort() for sorting compacted labels in %s at line %d\n", __FILE__, __LINE__);
			exit(-1);
		}

		//compute sigma
		res = cudppSegmentedScan(scanplan,
				d_sigsegment,
				d_dpixmi,
				d_scanmask,
				histo_size);
		if (CUDPP_SUCCESS != res)
		{
			printf("Error in cudppSort() for sorting compacted labels in %s at line %d\n", __FILE__, __LINE__);
			exit(-1);
		}

		// Destroy the segment prefix scan plan
		res = cudppDestroyPlan(scanplan);
		if (CUDPP_SUCCESS != res)
		{
			printf("Error in cudppSort() for sorting compacted labels in %s at line %d\n", __FILE__, __LINE__);
			exit(-1);
		}

		//////////////////////////////////////
		//compute the absdiff(lowsum and highsum)




		//update histo_size
		//update histo segment (compact)
		//update startpos
		n--;
	}

//	int grid1 = 4096;
//	int block1 = 1024;
//
//	backguess_kernel1<<<grid1, block1, 48*1024>>>(
//			d_nlevels,
//			d_mean,
//			d_sigma,
//			d_qscale,
//			d_qzero,
//			d_histo,
//			d_lowsum,
//			d_highsum,
//			nx_mesh*ny_mesh,
//			histo_size);

	/*
	int *h_lowsum = (int*)malloc(histo_size*sizeof(int));
	int *h_highsum  = (int*)malloc(histo_size*sizeof(int));
	int *h_level = (int*)malloc(nx_mesh*ny_mesh*sizeof(int));

	cudaMemcpy(h_lowsum, d_lowsum, histo_size*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_highsum, d_highsum, histo_size*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_level, d_nlevels, nx_mesh*ny_mesh*sizeof(int), cudaMemcpyDeviceToHost);

	int mindiff = INT_MAX;
	int minpos;
	for(int i=0; i<h_level[0]; i++)
	{
		if(abs(h_lowsum[i]-h_highsum[i]) < mindiff)
		{
			mindiff = abs(h_lowsum[i]-h_highsum[i]);
			minpos = i;
		}

		//if(i < 4000)
		//	printf("lowsum , highsum is %d\t%d\n", h_lowsum[i], h_highsum[i]);
	}

	printf("pos and lowsum, highsum is %d, %d, %d \n", minpos, h_lowsum[minpos], h_highsum[minpos]);

	free(h_level);
	free(h_lowsum);
	free(h_highsum);
	*/

	cudaFree(d_scanmask);
	cudaFree(d_lowsum);
	cudaFree(d_highsum);
	cudaFree(d_histosegment);

	cudaFree(d_dpix);
	cudaFree(d_dpixmi);
	cudaFree(d_meansegment);
	cudaFree(d_sigsegment);

	cudaFree(d_startpos);
	cudaFree(d_lcut);
	cudaFree(d_hcut);
	cudaFree(d_sig);
	cudaFree(d_sig1);

}

/**
 * "Look for `bad' meshes and interpolate them if necessary."
 * This special case part is not implemented currently.
 */
__global__ void medianfilter_kernel(
		float *d_mean,
		float *d_sigma,
		float *d_outmean,
		float *d_outsigma,
		double backfthresh,
		int npx,
		int npy,
		int nx_mesh,
		int ny_mesh)
{
	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	if(tid_x >= nx_mesh || tid_y >= ny_mesh)
		return;

	const int tid = tid_y * nx_mesh + tid_x;

	int mask_width, mask_height, mask_size;

	//read current data and its neighbours
	if(tid_x <= npx)
		mask_width = tid_x*2 + 1;
	else if(tid_x >= nx_mesh - npx)
		mask_width = (nx_mesh - tid_x - 1)*2 + 1;
	else
		mask_width = npx*2 + 1;

	if(tid_y <= npy)
		mask_height = tid_y*2 + 1;
	else if(tid_y >= ny_mesh - npy)
		mask_height = (ny_mesh - tid_y -1)*2 + 1;
	else
		mask_height = npy*2 + 1;

	mask_size = mask_width * mask_height;

	float *mean_array = (float*)malloc(mask_size*sizeof(float));
	float *sigma_array = (float*)malloc(mask_size*sizeof(float));

	int k = 0;

	for(int i=tid_x - mask_width/2; i <=tid_x + mask_width/2; i++)
		for(int j=tid_y - mask_height/2; j<=tid_y + mask_height/2; j++)
		{
			mean_array[k] = d_mean[j*nx_mesh+i];
			sigma_array[k++] = d_sigma[j*nx_mesh+i];
		}

	float median;
	if (fabs((median=hmedian(mean_array, mask_size))-d_mean[tid])>=backfthresh)
	{
		d_outmean[tid] = median;
		d_outsigma[tid] = hmedian(sigma_array, mask_size);
	}
	else
	{
		d_outmean[tid] = d_mean[tid];
		d_outsigma[tid] = d_sigma[tid];
	}
}

/**
 * Pre-compute 2nd derivatives along the y direction at background nodes.
 * shape: one dimension, one thread handle a column
 *
 * d_in:  back or sigma
 * d_out: dback or dsigma
 */
__global__ void makebackspline_kernel(
		float *d_mapin,
		float *d_mapout,
		int nx_mesh,
		int ny_mesh)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= nx_mesh)
		return;

	int nbym1 = ny_mesh - 1;

	float *mapt = d_mapin + tid;
	float *dmapt = d_mapout + tid;
	float *u, temp;

	if (ny_mesh > 1)
	{
		u = (float*)malloc(nbym1*sizeof(float)); /* temporary array */
		*dmapt = *u = 0.0;	/* "natural" lower boundary condition */
		mapt += nx_mesh;
		for (int y=1; y<nbym1; y++, mapt+=nx_mesh)
		{
			temp = -1/(*dmapt+4);
			*(dmapt += nx_mesh) = temp;
			temp *= *(u++) - 6*(*(mapt+nx_mesh)+*(mapt-nx_mesh)-2**mapt);
			*u = temp;
		}
		*(dmapt+=nx_mesh) = 0.0;	/* "natural" upper boundary condition */

		for (int y=ny_mesh-2; y--;)
		{
			temp = *dmapt;
			dmapt -= nx_mesh;
			*dmapt = (*dmapt*temp+*(u--))/6.0;
		}
		free(u);
	}
	else
		*dmapt = 0.0;
}

__global__ void subconstback_kernel(
		float *d_imgbuf,
		float backmean,
		int width,
		int height)
{
	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	if(tid_x >= width || tid_y >= height)
		return;

	d_imgbuf[tid_y*width+tid_x] -= backmean;
}

__global__ void subbackground_kernel(
		float *d_imgbuf,
		float *d_back,
		float *d_dback,
		float *node_t,
		float *dnode_t,
		float *u_t,
		int nx_mesh,
		int ny_mesh,
		int backsize,
		int width,
		int height)
{
	int		i,j,x,yl;
	float	dx,dx0,dy, cdx,cdy, temp, *blo,*bhi,*dblo,*dbhi ;

	float *nodep;

	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	float 	*line = d_imgbuf + tid*width;

	float *node = node_t + tid*nx_mesh;
	float *dnode = dnode_t + tid*nx_mesh;
	float *u = u_t + tid*nx_mesh;

	if(tid >= height)
		return;

	if (ny_mesh > 1)
	{
		dy = (float)tid/backsize - 0.5;
		dy -= (yl = (int)dy);

		if (yl<0)
		{
			yl = 0;
			dy -= 1.0;
		}
		else if (yl>=ny_mesh-1)
		{
			yl = ny_mesh<2 ? 0 : ny_mesh-2;
			dy += 1.0;
		}

		/*-- Interpolation along y for each node */
		cdy = 1 - dy;
		blo = d_back + nx_mesh*yl;
		bhi = blo + nx_mesh;
		dblo = d_dback + nx_mesh*yl;
		dbhi = dblo + nx_mesh;

		//node = (float*)malloc(nx_mesh*sizeof(float));
		nodep = node;
		for (x=nx_mesh; x--;)
			*(nodep++) = (cdy**(blo++) + dy**(bhi++) + (cdy*cdy*cdy-cdy)**(dblo++) + (dy*dy*dy-dy)**(dbhi++));

		/*-- Computation of 2nd derivatives along x */
		//dnode = (float*)malloc(nx_mesh*sizeof(float));	/* 2nd derivative along x */
		if (nx_mesh>1)
		{
			//u = (float*)malloc((nx_mesh-1)*sizeof(float)); /* temporary array */
			*dnode = *u = 0.0;	/* "natural" lower boundary condition */
			nodep = node+1;
			for (x=(nx_mesh-1); --x; nodep++)
			{
				temp = -1/(*(dnode++)+4);
				*dnode = temp;
				temp *= *(u++) - 6*(*(nodep+1)+*(nodep-1)-2**nodep);
				*u = temp;
			}

			*(++dnode) = 0.0;	/* "natural" upper boundary condition */
			for (x=nx_mesh-2; x--;)
			{
				temp = *(dnode--);
				*dnode = (*dnode*temp+*(u--))/6.0;
			}
			//free(u);
			dnode--;
		}
	}
	else
	{
		/*-- No interpolation and no new 2nd derivatives needed along y */
		node = d_back;
		dnode = d_dback;
	}

	/*-- Interpolation along x */
	if (nx_mesh > 1)
	{
		dx  = (1.0/backsize - 1)/2;	/* dx of the first pixel in the row */
		dx0 = ((backsize+1)%2)*(1.0/backsize)/2;	/* dx of the 1st pixel right to a bkgnd node */
		blo = node;
		bhi = node + 1;
		dblo = dnode;
		dbhi = dnode + 1;
		for (x=i=0,j=width; j--; i++, dx += (1.0/backsize))
		{
			if (i==backsize/2 && x>0 && x<(nx_mesh-1))
			{
				blo++;
				bhi++;
				dblo++;
				dbhi++;
				dx = dx0;
			}
			cdx = 1 - dx;
			*(line++) -= (float)(cdx*(*blo+(cdx*cdx-1)**dblo) + dx*(*bhi+(dx*dx-1)**dbhi));
			if (i==backsize)
			{
				x++;
				i = 0;
			}
		}
	}
	else
		for (j=width; j--;)
			*(line++) -= (float)*node;

	if (ny_mesh > 1)
	{
		//free(node);
		//free(dnode);
	}
}

__global__ void computeSqure(float *d_in, float *d_out, float *d_grid, int width, int height)
{
	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	if(tid_x < width && tid_y < height) {

		const int tid = tid_y * width + tid_x;
		d_out[tid] = d_in[tid] * d_in[tid];

		int x = tid_x/64;
		int y = tid_y/64;

		d_grid[tid] = (y*64+x);
	}
}

/**
 * only square back shape(backw=backh) is supported.
 */
extern "C" void back_estimate(
		int back_size,
		float backfthresh,
		float *backmean,
		float *backsigma,
		int nbackfx,
		int nbackfy,
		int back_type )
{

	if(back_size < SQUARE_BLK_WIDTH)
	{
		printf("bkgnd mesh size should be larger than block size.\n"
				"Try a smaller blk size for bkgnd estimation\n");
		exit(-1);
	}
	int nx_mesh = (width-1)/back_size+1;
	int ny_mesh = (height-1)/back_size+1;

	int 	*d_npix;
	int 	*d_nlevels;
	float 	*d_lcut;
	float 	*d_hcut;
	float	*d_qscale;
	float	*d_qzero;
	int		*d_histo;

	float	*d_mean1;
	float	*d_sigma1;

	float time;

	checkCudaErrors(cudaMalloc((void**)&d_npix, 	nx_mesh * ny_mesh * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&d_nlevels, 	nx_mesh * ny_mesh * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&d_mean, 	nx_mesh * ny_mesh * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_sigma, 	nx_mesh * ny_mesh * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_lcut, 	nx_mesh * ny_mesh * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_hcut, 	nx_mesh * ny_mesh * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_qscale, 	nx_mesh * ny_mesh * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_qzero, 	nx_mesh * ny_mesh * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_histo, 	nx_mesh * ny_mesh * QUANTIF_NMAXLEVELS*sizeof(int)));

	//store the median filter result
	checkCudaErrors(cudaMalloc((void**)&d_mean1, 	nx_mesh * ny_mesh * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_sigma1, 	nx_mesh * ny_mesh * sizeof(float)));

	cudaEventRecord(start, 0);

	dim3 grid(nx_mesh, ny_mesh);
	dim3 block(SQUARE_BLK_WIDTH, SQUARE_BLK_WIDTH);

	//48*1024 should be changed (21ms)
	int smem_size = back_size*back_size*sizeof(float) +
			(SQUARE_BLK_WIDTH*SQUARE_BLK_WIDTH)*(2*sizeof(double)+2*sizeof(int)) + 2*sizeof(float);
	int smem_size1 = back_size*back_size*sizeof(float) +
			QUANTIF_NMAXLEVELS * sizeof(int);

	smem_size = max(smem_size, smem_size1);

	if(smem_size > 48*1024)
	{
		printf("Two much share memory space needed, exit!");
		exit(-1);
	}

	//25 ms
	backstat_kernel<<<grid, block, smem_size>>>(d_pixelArray, d_npix, d_nlevels, d_mean, d_sigma,
			d_lcut, d_hcut, d_qscale, d_qzero, d_histo, width, height, back_size);


	dim3 grid1((nx_mesh-1)/SQUARE_BLK_WIDTH+1, (ny_mesh-1)/SQUARE_BLK_WIDTH+1);
	dim3 block1(SQUARE_BLK_WIDTH, SQUARE_BLK_WIDTH);

	//102~108 ms
	backguess_kernel<<<grid1, block1>>>(d_nlevels, d_mean, d_sigma, d_qscale,d_qzero, d_histo,
			/*d_mean1,d_sigma1,*/ nx_mesh, ny_mesh);

	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
	{
		printf("cudaerror %d\n", err);
	}

	int npx = nbackfx/2;
	int npy = nbackfy/2;

	/* Median-filter and check suitability of the background map */
	//NFPRINTF(OUTPUT, "Filtering background map(s)");

	//2 ms
	//we use the d_mean1 and d_sigma1 as the output array
	medianfilter_kernel<<<grid1, block1>>>(d_mean,
			d_sigma,
			d_mean1,
			d_sigma1,
			backfthresh,
			npx, npy,
			nx_mesh,
			ny_mesh);


	float *d_dback = d_lcut; //d_lcut is re-used to store d_dback
	float *d_dsigma = d_hcut; //d_hcut is re-used to store d_dsigma

	int grid2 = (nx_mesh-1)/MAX_THREADS_PER_BLK + 1;
	int block2 = MAX_THREADS_PER_BLK;

	/* Compute 2nd derivatives along the y-direction */
	//NFPRINTF(OUTPUT, "Computing background d-map"); 2 ms
	makebackspline_kernel<<<grid2, block2>>>(d_mean1, d_dback, nx_mesh, ny_mesh);

	//NFPRINTF(OUTPUT, "Computing background-noise d-map"); 2 ms
	makebackspline_kernel<<<grid2, block2>>>(d_sigma1,d_dsigma,nx_mesh, ny_mesh);


	//here we use d_mean and d_sigma to store the backup mean and sigma
	//since the element order will be changed after computing median mean and sigma
	checkCudaErrors(cudaMemcpy(d_mean,  d_mean1, nx_mesh*ny_mesh*sizeof(float),
			cudaMemcpyDeviceToDevice));
	checkCudaErrors(cudaMemcpy(d_sigma, d_sigma1, nx_mesh*ny_mesh*sizeof(float),
			cudaMemcpyDeviceToDevice));

	//0.57 ms
	//this method will change the content(order) of d_mean and d_sigma array
	get_mean_sigma(d_mean1, d_sigma1, backmean, backsigma, nx_mesh*ny_mesh);

	/*-- In absolute background mode, just subtract a cste */
	if(back_type == 1)
	{
		dim3 grid3((width - 1) / SQUARE_BLK_WIDTH + 1,
				(height - 1) / SQUARE_BLK_HEIGHT + 1);
		dim3 block3(SQUARE_BLK_WIDTH, SQUARE_BLK_HEIGHT);

		subconstback_kernel<<<grid3, block3>>>(d_pixelArray, *backmean, width, height);
	}else
	{
		int grid4 = (height-1)/MAX_THREADS_PER_BLK + 1;
		int block4= MAX_THREADS_PER_BLK;

		float *node_t, *dnode_t, *u_t;
		cudaMalloc((void**)&node_t,  height*nx_mesh*sizeof(float));
		cudaMalloc((void**)&dnode_t, height*nx_mesh*sizeof(float));
		cudaMalloc((void**)&u_t, height*nx_mesh*sizeof(float));

		//137 ms (35 ms if remove dynamic allocation parts)
		/* this kernel function may have precision problem */
		subbackground_kernel<<<grid4, block4>>>(d_pixelArray,
				d_mean,
				d_dback,
				node_t,
				dnode_t,
				u_t,
				nx_mesh,
				ny_mesh,
				back_size,
				width,
				height);

		cudaFree(node_t);
		cudaFree(dnode_t);
		cudaFree(u_t);
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

#ifdef DEBUG_CUDA
	printf("time consumed by cuda back estimation %f\n", time);
#endif

	cudaFree(d_npix);
	cudaFree(d_nlevels);
	cudaFree(d_lcut);
	cudaFree(d_hcut);
	cudaFree(d_qscale);
	cudaFree(d_qzero);
	cudaFree(d_histo);
	cudaFree(d_mean1);
	cudaFree(d_sigma1);

}

void get_mean_sigma(float *d_mean1,
		float *d_sigma1,
		float *backmean,
		float *backsig,
		int nb)
{
	CUDPPResult res;

	config.datatype = CUDPP_FLOAT;
	config.algorithm = CUDPP_SORT_RADIX;
	config.options = CUDPP_OPTION_KEYS_ONLY;

	//create the cudpp sort plan
	res = cudppPlan(theCudpp, &sortplan, config, nb, 1, 0);
	if (CUDPP_SUCCESS != res)
	{
		printf("Error creating CUDPP sort Plan in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	res = cudppRadixSort(sortplan,  d_mean1, NULL, nb);
	if (CUDPP_SUCCESS != res)
	{
		printf("Error in cudppSort() for sorting mean array in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	res = cudppRadixSort(sortplan,  d_sigma1, NULL, nb);
	if (CUDPP_SUCCESS != res)
	{
		printf("Error in cudppSort() for sorting mean array in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	// Destroy the sort plan
	res = cudppDestroyPlan(sortplan);
	if (CUDPP_SUCCESS != res)
	{
		printf("Error destroying CUDPP sort Plan in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	if(nb%2 == 0)
	{
		float median[2];
		float sigma[2];
		cudaMemcpy(median,(d_mean1 + nb/2-1),  2*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(sigma, (d_sigma1 + nb/2-1), 2*sizeof(float), cudaMemcpyDeviceToHost);

		*backmean = (median[0]+median[1])/2;
		*backsig = (sigma[0]+sigma[1])/2;
	}
	else
	{
		cudaMemcpy(backmean,(d_mean1 + nb/2),  sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(backsig,(d_sigma1 + nb/2),   sizeof(float), cudaMemcpyDeviceToHost);
	}
}

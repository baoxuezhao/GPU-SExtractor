/*
 * cudainit.cu
 *
 *  Created on: 10 Jan, 2013
 *      Author: zhao
 */
#include <stdio.h>
#include <cuda.h>
#include <time.h>

#include "../util5/cudpp.h"
#include "../util5/helper_cuda.h"
#include "cudatypes.h"


extern "C" void init_device(int _width, int _height, float *imgbuf)
{
	//cudaDeviceReset();

	width = _width;
	height = _height;

	cudaEventCreate(&start_t);
	cudaEventCreate(&stop_t);

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start_t, 0);


	checkCudaErrors(cudaMalloc((void**)&d_cdPixArray, width*height * sizeof(float)));

	checkCudaErrors(cudaMalloc((void**)&d_pixelArray, width*height * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&d_labelArray, width*height * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&d_equivArray, width*height * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&d_compactMask,width*height * sizeof(unsigned int)));

	checkCudaErrors(cudaMalloc((void**)&d_numPixAboveThresh, sizeof(size_t)));


	float time;
	cudaEventRecord(start, 0);

	cudaMemcpy((float*)d_pixelArray, imgbuf, width*height*sizeof(float), cudaMemcpyHostToDevice);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

#ifdef DEBUG_CUDA
	printf("time consumed by cuda data transfer in is %f\n", time);
#endif

	//Initialize the CUDPP Library
	cudppCreate(&theCudpp);

}

void init_objects(int total)
{
	d_index 	= d_cuttedObjLabelArray;
	d_flag 		= d_cuttedObjFlagArray;
	d_fdnpix 	= d_cuttedPixCountArray;
	d_dthresh 	= d_cuttedDthreshArray;

	cudaMalloc((void**)&d_xmin, total*sizeof(unsigned int));
	cudaMalloc((void**)&d_xmax, total*sizeof(unsigned int));
	cudaMalloc((void**)&d_ymin, total*sizeof(unsigned int));
	cudaMalloc((void**)&d_ymax, total*sizeof(unsigned int));
	cudaMalloc((void**)&d_dnpix,total*sizeof(unsigned int));
	cudaMalloc((void**)&d_npix, total*sizeof(unsigned int));
	cudaMalloc((void**)&d_peakx,total*sizeof(unsigned int));
	cudaMalloc((void**)&d_peaky,total*sizeof(unsigned int));

	cudaMalloc((void**)&d_bkg, 		total*sizeof(float));
	cudaMalloc((void**)&d_dbkg, 	total*sizeof(float));
	cudaMalloc((void**)&d_sigbkg, 	total*sizeof(float));
	cudaMalloc((void**)&d_a, 		total*sizeof(float));
	cudaMalloc((void**)&d_b, 		total*sizeof(float));
	cudaMalloc((void**)&d_cxx,		total*sizeof(float));
	cudaMalloc((void**)&d_cxy,		total*sizeof(float));
	cudaMalloc((void**)&d_cyy,		total*sizeof(float));
	cudaMalloc((void**)&d_theta,	total*sizeof(float));
	cudaMalloc((void**)&d_abcor,	total*sizeof(float));
	cudaMalloc((void**)&d_peak, 	total*sizeof(float));
	cudaMalloc((void**)&d_dpeak, 	total*sizeof(float));
	cudaMalloc((void**)&d_fdpeak, 	total*sizeof(float));
	cudaMalloc((void**)&d_flux, 	total*sizeof(float));
	cudaMalloc((void**)&d_dflux,	total*sizeof(float));
	cudaMalloc((void**)&d_fdflux,	total*sizeof(float));
	cudaMalloc((void**)&d_fluxerr,	total*sizeof(float));
	cudaMalloc((void**)&d_thresh,	total*sizeof(float));
	cudaMalloc((void**)&d_mthresh,	total*sizeof(float));
	cudaMalloc((void**)&d_fwhm, 	total*sizeof(float));

	cudaMalloc((void**)&d_mx, 		total*sizeof(double));
	cudaMalloc((void**)&d_my, 		total*sizeof(double));
	cudaMalloc((void**)&d_mx2, 		total*sizeof(double));
	cudaMalloc((void**)&d_my2, 		total*sizeof(double));
	cudaMalloc((void**)&d_mxy, 		total*sizeof(double));
	cudaMalloc((void**)&d_poserr_mx2, 	total*sizeof(double));
	cudaMalloc((void**)&d_poserr_my2, 	total*sizeof(double));
	cudaMalloc((void**)&d_poserr_mxy, 	total*sizeof(double));
	cudaMalloc((void**)&d_singuflag,	total*sizeof(char));

	//for astrom measurement
	cudaMalloc((void**)&d_mxw, 		total*sizeof(double));
	cudaMalloc((void**)&d_myw, 		total*sizeof(double));
	cudaMalloc((void**)&d_alphas, 	total*sizeof(double));
	cudaMalloc((void**)&d_deltas, 	total*sizeof(double));
	cudaMalloc((void**)&d_alpha2000,total*sizeof(double));
	cudaMalloc((void**)&d_delta2000,total*sizeof(double));

	//cudaMalloc((void**)&d_iso,			NISO*sizeof(int*));
	for(int i=0; i<NISO; i++)
	{
		cudaMalloc((void**)&d_iso[i],	total*sizeof(int));
		cudaMemset(d_iso[i], 0, total*sizeof(int));
	}

	//for variables used in endobject phase
	cudaMalloc((void**)&d_posx, total*sizeof(double));
	cudaMalloc((void**)&d_posy, total*sizeof(double));

	cudaMalloc((void**)&d_elong, 	total*sizeof(float));
	cudaMalloc((void**)&d_ellip, 	total*sizeof(float));
	cudaMalloc((void**)&d_polar, 	total*sizeof(float));
	cudaMalloc((void**)&d_sprob, 	total*sizeof(float));
	cudaMalloc((void**)&d_sposx, 	total*sizeof(float));
	cudaMalloc((void**)&d_sposy, 	total*sizeof(float));
	//cudaMalloc((void**)&d_poserr_a, total*sizeof(float));
	//cudaMalloc((void**)&d_poserr_b, total*sizeof(float));
	//cudaMalloc((void**)&d_poserr_theta, total*sizeof(float));
	//cudaMalloc((void**)&d_poserr_cxx, total*sizeof(float));
	//cudaMalloc((void**)&d_poserr_cyy, total*sizeof(float));
	//cudaMalloc((void**)&d_poserr_cxy, total*sizeof(float));
	cudaMalloc((void**)&d_flux_iso, total*sizeof(float));
	cudaMalloc((void**)&d_fluxerr_iso, total*sizeof(float));
	cudaMalloc((void**)&d_flux_isocor, total*sizeof(float));
	cudaMalloc((void**)&d_fluxerr_isocor, total*sizeof(float));
	cudaMalloc((void**)&d_mag_iso, total*sizeof(float));
	cudaMalloc((void**)&d_magerr_iso, total*sizeof(float));
	cudaMalloc((void**)&d_kronfactor, total*sizeof(float));
	cudaMalloc((void**)&d_flux_auto, total*sizeof(float));
	cudaMalloc((void**)&d_fluxerr_auto, total*sizeof(float));
	cudaMalloc((void**)&d_mag_auto, total*sizeof(float));
	cudaMalloc((void**)&d_magerr_auto, total*sizeof(float));

	for(int i=0; i<NAPER; i++)
	{
		cudaMalloc((void**)&d_flux_aper[i],	total*sizeof(float));
		cudaMalloc((void**)&d_fluxerr_aper[i],	total*sizeof(float));
		cudaMalloc((void**)&d_mag_aper[i],	total*sizeof(float));
		cudaMalloc((void**)&d_magerr_aper[i],	total*sizeof(float));
	}

}


extern "C" void clear_device(float *img)
{
	cudaMemcpy(img, d_pixelArray, width*height*sizeof(float), cudaMemcpyDeviceToHost);

	checkCudaErrors(cudaFree(d_cdPixArray));
	checkCudaErrors(cudaFree(d_pixelArray));

	checkCudaErrors(cudaFree(d_labelArray));
	checkCudaErrors(cudaFree(d_equivArray));
	checkCudaErrors(cudaFree(d_compactMask));
	checkCudaErrors(cudaFree(d_numPixAboveThresh));

	checkCudaErrors(cudaFree(d_compactedIndexArray));
	checkCudaErrors(cudaFree(d_compactedParentLabel));
	//checkCudaErrors(cudaFree(d_compactedPixelArray));
	checkCudaErrors(cudaFree(d_compactedcdPixelArray));
	checkCudaErrors(cudaFree(d_segmentMask));

	checkCudaErrors(cudaFree(d_pixelCountMask));
	checkCudaErrors(cudaFree(d_pixelCountSegment));
	checkCudaErrors(cudaFree(d_fdpeakSegment));
	checkCudaErrors(cudaFree(d_fdfluxSegment));
	checkCudaErrors(cudaFree(d_prunedSegmentMask));

	/*-----------data structure to store objects after cutting-----------*/
	//unfinished, need to free all
	checkCudaErrors(cudaFree(d_cuttedObjLevelArray));
	checkCudaErrors(cudaFree(d_cuttedObjLabelArray));
	checkCudaErrors(cudaFree(d_cuttedObjIndexArray));
	checkCudaErrors(cudaFree(d_cuttedObjFlagArray));
	checkCudaErrors(cudaFree(d_cuttedPixCountArray));
	checkCudaErrors(cudaFree(d_cuttedRootlabelArray));
	checkCudaErrors(cudaFree(d_cuttedDthreshArray));

	checkCudaErrors(cudaFree(d_xmin));
	checkCudaErrors(cudaFree(d_xmax));
	checkCudaErrors(cudaFree(d_ymin));
	checkCudaErrors(cudaFree(d_ymax));
	checkCudaErrors(cudaFree(d_dnpix));
	checkCudaErrors(cudaFree(d_mx));
	checkCudaErrors(cudaFree(d_my));
	checkCudaErrors(cudaFree(d_mx2));
	checkCudaErrors(cudaFree(d_my2));
	checkCudaErrors(cudaFree(d_mxy));
	checkCudaErrors(cudaFree(d_cxx));
	checkCudaErrors(cudaFree(d_cxy));
	checkCudaErrors(cudaFree(d_cyy));
	checkCudaErrors(cudaFree(d_a));
	checkCudaErrors(cudaFree(d_b));
	checkCudaErrors(cudaFree(d_theta));
	checkCudaErrors(cudaFree(d_abcor));
	checkCudaErrors(cudaFree(d_fdpeak));
	checkCudaErrors(cudaFree(d_dpeak));
	checkCudaErrors(cudaFree(d_fdflux));
	checkCudaErrors(cudaFree(d_dflux));
	checkCudaErrors(cudaFree(d_singuflag));

	for(int i=0; i<NISO; i++)
		checkCudaErrors(cudaFree(d_iso[i]));

	checkCudaErrors(cudaFree(d_posx));
	checkCudaErrors(cudaFree(d_posy));
	checkCudaErrors(cudaFree(d_elong));
	checkCudaErrors(cudaFree(d_ellip));
	checkCudaErrors(cudaFree(d_polar));
	checkCudaErrors(cudaFree(d_sprob));
	checkCudaErrors(cudaFree(d_sposx));
	checkCudaErrors(cudaFree(d_sposy));
	checkCudaErrors(cudaFree(d_flux_iso));
	checkCudaErrors(cudaFree(d_fluxerr_iso));
	checkCudaErrors(cudaFree(d_flux_isocor));
	checkCudaErrors(cudaFree(d_fluxerr_isocor));
	checkCudaErrors(cudaFree(d_mag_iso));
	checkCudaErrors(cudaFree(d_magerr_iso));
	checkCudaErrors(cudaFree(d_kronfactor));
	checkCudaErrors(cudaFree(d_flux_auto));
	checkCudaErrors(cudaFree(d_fluxerr_auto));
	checkCudaErrors(cudaFree(d_mag_auto));
	checkCudaErrors(cudaFree(d_magerr_auto));

	for(int i=0; i<NAPER; i++)
	{
		checkCudaErrors(cudaFree(d_flux_aper[i]));
		checkCudaErrors(cudaFree(d_fluxerr_aper[i]));
		checkCudaErrors(cudaFree(d_mag_aper[i]));
		checkCudaErrors(cudaFree(d_magerr_aper[i]));
	}

	//for astrom measurement
	checkCudaErrors(cudaFree(d_mxw));
	checkCudaErrors(cudaFree(d_myw));
	checkCudaErrors(cudaFree(d_alphas));
	checkCudaErrors(cudaFree(d_deltas));
	checkCudaErrors(cudaFree(d_alpha2000));
	checkCudaErrors(cudaFree(d_delta2000));

	/*---------------data structure to background result-----------*/
	checkCudaErrors(cudaFree(d_mean));
	checkCudaErrors(cudaFree(d_sigma));

	checkCudaErrors(cudaFree(d_masterIndex));

	// Shut down the CUDPP library
	cudppDestroy(theCudpp);

	/*-----------------timing variables-------------------*/
	cudaEventRecord(stop_t, 0);
	cudaEventSynchronize(stop_t);

	float totaltime;
	cudaEventElapsedTime(&totaltime, start_t, stop_t);

	cudaEventDestroy(start_t);
	cudaEventDestroy(stop_t);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	printf("Time counsumed by gpu code in SExtractor is: %f (ms)\n", totaltime);

}




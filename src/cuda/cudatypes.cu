/*
 * cudavars.cu
 *
 *  Created on: 10 Jan, 2013
 *      Author: zhao
 */
#include "../util5/cudpp.h"
#include "../util5/helper_cuda.h"
#include "cudatypes.h"

//corresponds to the cdscan field in sextractor
float			*d_cdPixArray; //cdvalue field

//corresponds to the the scan field in sextractor
float			*d_pixelArray; //dvalue field
int				*d_labelArray;
int				*d_equivArray;
unsigned int	*d_indexArray;

////////////////used to compact before detection
unsigned int	*d_compactMask;
size_t			*d_numPixAboveThresh;

////////////////temp variables, allocated in the first detection pass(except. d_compactedParentLabel)
unsigned int	*d_segmentMask;
unsigned int	*d_compactedIndexArray;
unsigned int	*d_compactedLabelArray;
unsigned int	*d_compactedParentLabel;
float			*d_compactedcdPixelArray;
float			*d_compactedPixelArray;

unsigned int	*d_pixelCountMask;
unsigned int	*d_pixelCountSegment;

float			*d_fdpeakSegment;
float			*d_fdfluxSegment;
unsigned int 	*d_prunedSegmentMask;

/*---------------The output of the first scan-------------*/
unsigned int	*d_pixelCountArray; //fdnpix
float			*d_fdpeakArray;
float			*d_fdfluxArray;
float			*d_dthreshArray;
unsigned int	*d_ok;

unsigned int	*d_finalPixelIndexArray;
unsigned int	*d_finalLabelArray;
unsigned int	*d_finalObjIndexArray;

/*---------------data structure to background result-----------*/
float 			*d_mean;
float 			*d_sigma;


/*-----------data structure to store objects after cutting-----------*/
unsigned int	*d_cuttedObjLevelArray;
unsigned int	*d_cuttedObjLabelArray;
unsigned int	*d_cuttedObjIndexArray;
unsigned int	*d_cuttedObjFlagArray;
unsigned int	*d_cuttedPixCountArray;
unsigned int	*d_cuttedRootlabelArray;
float			*d_cuttedDthreshArray;

/*-----------data structure to store objects attributes-----------*/
unsigned int	*d_index;
unsigned int	*d_flag;
unsigned int	*d_fdnpix;
float			*d_dthresh;

unsigned int 	*d_xmin;
unsigned int 	*d_xmax;
unsigned int 	*d_ymin;
unsigned int 	*d_ymax;
unsigned int 	*d_dnpix;
unsigned int 	*d_npix;
unsigned int	*d_peakx;
unsigned int	*d_peaky;

float 		*d_bkg;
float 		*d_dbkg;
float		*d_sigbkg;
float 		*d_a;
float 		*d_b;
float 		*d_cxx;
float 		*d_cxy;
float 		*d_cyy;
float 		*d_theta;
float 		*d_abcor;
float		*d_peak;
float 		*d_dpeak;
float 		*d_fdpeak;
float		*d_flux;
float 		*d_dflux;
float 		*d_fdflux;
float		*d_fluxerr;
float		*d_thresh;
float		*d_mthresh;
float		*d_fwhm;

double 		*d_mx;
double 		*d_my;
double 		*d_mx2;
double 		*d_my2;
double 		*d_mxy;
double		*d_poserr_mx2;
double		*d_poserr_my2;
double		*d_poserr_mxy;
char  		*d_singuflag;

int			*d_iso[NISO];
//////////////////////////
//used in endobject phase
double		*d_mxw;
double		*d_myw;
double		*d_alphas;
double		*d_deltas;
double		*d_alpha2000;
double		*d_delta2000;

double 		*d_posx;
double		*d_posy;
float		*d_elong;
float		*d_ellip;
float		*d_polar;
float		*d_sprob;
float		*d_sposx;
float		*d_sposy;
//float		*d_poserr_a;
//float		*d_poserr_b;
//float		*d_poserr_theta;
//float		*d_poserr_cxx;
//float		*d_poserr_cyy;
//float		*d_poserr_cxy;
float		*d_flux_iso;
float		*d_fluxerr_iso;
float		*d_flux_isocor;
float		*d_fluxerr_isocor;
float		*d_mag_iso;
float		*d_magerr_iso;
float		*d_kronfactor;
float		*d_flux_auto;
float		*d_fluxerr_auto;
float		*d_mag_auto;
float		*d_magerr_auto;

float		*d_flux_aper[NAPER];
float		*d_fluxerr_aper[NAPER];
float		*d_mag_aper[NAPER];
float		*d_magerr_aper[NAPER];

/*-----------CUDPP variables-----------*/
CUDPPHandle 	theCudpp;

CUDPPConfiguration config;
CUDPPHandle 	scanplan;
CUDPPHandle		reduceplan;
CUDPPHandle		segscanplan;
CUDPPHandle		compactplan;
CUDPPHandle		sortplan;

int 	width;
int 	height;
float	thresh;
float 	global_dthresh;

/*-----------timing variables-----------*/
cudaEvent_t start_t, stop_t;	//total time measurement
cudaEvent_t start, stop;		//time measurement of each part

unsigned int 	*d_masterIndex;

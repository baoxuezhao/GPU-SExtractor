/*
 				scan.c

 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	Part of:	SExtractor
 *
 *	Author:		E.BERTIN (IAP)
 *
 *	Contents:	functions for extraction of connected pixels from
 *			a pixmap.
 *
 *	Last modify:	11/01/2008
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"back.h"
#include	"check.h"
#include	"clean.h"
#include	"extract.h"
#include	"filter.h"
#include	"image.h"
#include	"plist.h"

#include 	"cuda/cudadetection.h"
#include 	"cuda/cudadeblend.h"
#include	"cuda/cudafilter.h"
#include	"cuda/cudaclean.h"
#include	"cuda/cudaback.h"
#include	"cuda/cudainit.h"
#include	"cuda/cudaanalyse.h"

extern tabstruct	*objtab;
extern FILE			*ascfile;
void printobj(FILE	*stream, outobjstruct *outobj);

/****************************** scanimage ************************************
PROTO   void scanimage(picstruct *field, picstruct *dfield, picstruct *ffield,
        picstruct *wfield, picstruct *dwfield)
PURPOSE Scan of the large pixmap(s). Main loop and heart of the program.
INPUT   Measurement field pointer,
        Detection field pointer,
        Flag field pointer,
        Measurement weight-map field pointer,
        Detection weight-map field pointer,
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 11/01/2008
 ***/
void	scanimage(picstruct *field, picstruct **pffield, int nffield)

{
	int i, numdetected;
	size_t numValidObj;
	size_t numValidPix;

	//dthresh: 270.12
	//thresh: 900.40

	field->ymin = 0;

	//init_device(field->width, field->height);
	init_detection();

	//here(before convolution, the two field strips are the same

	if(thefilter->convw != thefilter->convh)
		printf("convolution mask is not square, not supported by the current kernel");

	cudaFilter(thefilter->conv, thefilter->convw);

	run_detection(field->dthresh, field->thresh, prefs.ext_minarea, &numValidObj, &numValidPix);

	init_deblend(field->dthresh, prefs.deblend_nthresh, numValidObj, numValidPix);
	numdetected = parcel_out(field->dthresh, prefs.deblend_nthresh, prefs.deb_maxarea,
			numValidObj, numValidPix, prefs.deblend_mincont);
	clear_deblend(prefs.deblend_nthresh);

	analyse_gpu(
			flagobj.poserr_mx2,
			flagobj.peakx,
			flagobj.iso[0],
			flagobj.fwhm,
			field->backsig,
			field->gain,
			field->ngamma,
			field->thresh,
			field->dthresh,
			field->satur_level,
			plistexist_var,
			plistoff_var,
			prefs.weightgain_flag,
			prefs.ext_minarea,
			prefs.clean_flag,
			prefs.detect_type,
			field->nbackx,
			field->nbacky,
			field->backw,
			field->backh,
			field->width,
			numdetected);

	//used as the flag for output
	unsigned int *masterIndex = (unsigned int*)malloc(numdetected*sizeof(unsigned int));

	run_clean(masterIndex, prefs.clean_param, prefs.clean_stacksize, numdetected);

	//object measurement
	outobjstruct *objlist = (outobjstruct*)malloc(numdetected*sizeof(outobjstruct));

	endobject_gpu(
			&prefs,
			brain,
			&flagobj2,
			objlist,
			field->wcs,
			/*---field vars---*/
			field->ymin,
			field->ymax,
			field->width,
			field->height,
			field->backsig,
			field->thresh,
			field->satur_level,
			field->pixscale,
			field->ngamma,
			field->gain,
			numdetected);

	clear_detection();
	clear_device(field->strip);

	int number = 0;
	//object output
	for(i = 0; i<numdetected; i++)
	{
		if(masterIndex[i]==-1)
		{
			(objlist+i)->number = number++;
			printobj(ascfile, objlist+i);
		}
	}

	printf("Objects: detected %d \t sextracted %d\n", numdetected, number);

	free(objlist);

	return;

	//code blow in this function is not used.
	/////////////////////////////////////////////////////////////////////////
	objstruct *objectlist = (objstruct*)malloc(numdetected*sizeof(objstruct));

	int total = numdetected;

	unsigned int *h_index = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_flag = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_fdnpix = (unsigned int*)malloc(total*sizeof(unsigned int));
	float 		 *h_dthresh = (float*)malloc(total*sizeof(float));

	unsigned int *h_xmin  = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_xmax  = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_ymin  = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_ymax  = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_dnpix = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_npix  = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_peakx = (unsigned int*)malloc(total*sizeof(unsigned int));
	unsigned int *h_peaky = (unsigned int*)malloc(total*sizeof(unsigned int));

	float 	*h_bkg  = (float*)malloc(total*sizeof(float));
	float 	*h_dbkg  = (float*)malloc(total*sizeof(float));
	float	*h_sigbkg  = (float*)malloc(total*sizeof(float));
	float 	*h_a  = (float*)malloc(total*sizeof(float));
	float 	*h_b  = (float*)malloc(total*sizeof(float));
	float 	*h_cxx  = (float*)malloc(total*sizeof(float));
	float 	*h_cxy  = (float*)malloc(total*sizeof(float));
	float 	*h_cyy  = (float*)malloc(total*sizeof(float));
	float 	*h_theta  = (float*)malloc(total*sizeof(float));
	float 	*h_abcor  = (float*)malloc(total*sizeof(float));
	float	*h_peak  = (float*)malloc(total*sizeof(float));
	float 	*h_dpeak  = (float*)malloc(total*sizeof(float));
	float 	*h_fdpeak  = (float*)malloc(total*sizeof(float));
	float	*h_flux  = (float*)malloc(total*sizeof(float));
	float 	*h_dflux = (float*)malloc(total*sizeof(float));
	float 	*h_fdflux  = (float*)malloc(total*sizeof(float));
	float	*h_fluxerr = (float*)malloc(total*sizeof(float));
	float	*h_thresh  = (float*)malloc(total*sizeof(float));
	float	*h_mthresh  = (float*)malloc(total*sizeof(float));
	float	*h_fwhm  = (float*)malloc(total*sizeof(float));

	double 	*h_mx   = (double*)malloc(total*sizeof(double));
	double 	*h_my   = (double*)malloc(total*sizeof(double));
	double 	*h_mx2   = (double*)malloc(total*sizeof(double));
	double 	*h_my2   = (double*)malloc(total*sizeof(double));
	double 	*h_mxy   = (double*)malloc(total*sizeof(double));
	double	*h_poserr_mx2   = (double*)malloc(total*sizeof(double));
	double	*h_poserr_my2   = (double*)malloc(total*sizeof(double));
	double	*h_poserr_mxy   = (double*)malloc(total*sizeof(double));
	char  	*h_singuflag   = (char*)malloc(total*sizeof(char));
	int		*h_iso[NISO];

	for(i=0; i<NISO; i++)
		h_iso[i] = (int*)malloc(total*sizeof(int));

	copyobjects(
			h_index,
			h_flag,
			h_fdnpix,
			h_dthresh,
			h_xmin,
			h_xmax,
			h_ymin ,
			h_ymax,
			h_dnpix,
			h_npix,
			h_peakx,
			h_peaky ,

			h_bkg,
			h_dbkg ,
			h_sigbkg,
			h_a,
			h_b,
			h_cxx,
			h_cxy,
			h_cyy,
			h_theta,
			h_abcor,
			h_peak,
			h_dpeak,
			h_fdpeak ,
			h_flux ,
			h_dflux,
			h_fdflux,
			h_fluxerr,
			h_thresh,
			h_mthresh,
			h_fwhm,
			h_mx,
			h_my,
			h_mx2,
			h_my2,
			h_mxy,
			h_poserr_mx2 ,
			h_poserr_my2,
			h_poserr_mxy,
			h_singuflag,
			h_iso,
			total);


	int idx, t;
	for(idx=0; idx<total; idx++)
	{
		objstruct* obj = objectlist+idx;

		obj->number = h_index[idx];
		obj->flag = h_flag[idx];
		obj->fdnpix = h_fdnpix[idx]; //v

		//if(idx < 10)
		//	printf("distance is %d\n", obj->fdnpix);

		obj->dthresh = h_dthresh[idx]; //v
		obj->xmin = h_xmin[idx]; //v
		obj->xmax = h_xmax[idx]; //v
		obj->ymin = h_ymin[idx]; //v
		obj->ymax = h_ymax[idx]; //v
		obj->dnpix = h_dnpix[idx]; //v
		obj->npix = h_npix[idx]; //v
		obj->peakx = h_peakx[idx]; //v
		obj->peaky = h_peaky[idx]; //v

		obj->bkg = h_bkg[idx]; //v
		obj->dbkg = h_dbkg[idx]; //v
		obj->sigbkg = h_sigbkg[idx]; //v
		obj->a = h_a[idx]; //v
		obj->b = h_b[idx]; //v
		obj->cxx = h_cxx[idx]; //v
		obj->cxy = h_cxy[idx]; //v
		obj->cyy = h_cyy[idx]; //v
		obj->theta = h_theta[idx];
		obj->abcor = h_abcor[idx];
		obj->peak = h_peak[idx];
		obj->dpeak = h_dpeak[idx];
		obj->fdpeak = h_fdpeak[idx];
		obj->flux = h_flux[idx];
		obj->dflux = h_dflux[idx];
		obj->fdflux = h_fdflux[idx];
		obj->fluxerr = h_fluxerr[idx];
		obj->thresh = h_thresh[idx];
		obj->mthresh = h_mthresh[idx];
		obj->fwhm = h_fwhm[idx];

		obj->mx =h_mx[idx];
		obj->my = h_my[idx];
		obj->mx2 = h_mx2[idx];
		obj->my2 = h_my2[idx];
		obj->mxy = h_mxy[idx];
		obj->poserr_mx2 = h_poserr_mx2[idx];
		obj->poserr_my2 = h_poserr_my2[idx];
		obj->poserr_mxy = h_poserr_mxy[idx];
		obj->singuflag = h_singuflag[idx];

		for(t=0; t<NISO; t++)
			obj->iso[t] = h_iso[t][idx];
	}

	printf("objectlist[1] %f, %d \n", objectlist[1].flux, objectlist[1].fdnpix);
	//mergeobject, pasteimage


	return;
}

void printobj(FILE	*stream, outobjstruct *outobj)
{
	fprintf(stream, "%10d", outobj->number);
	putc(' ', stream);
	fprintf(stream, "%11.7f", outobj->alpha2000);
	putc(' ', stream);
	fprintf(stream, "%+11.7f", outobj->delta2000);
	putc(' ', stream);

	fprintf(stream, "%10.3f",  outobj->sposx);
	putc(' ', stream);
	fprintf(stream, "%10.3f",  outobj->sposy);
	putc(' ', stream);
	fprintf(stream, "%12.7g",  outobj->flux_iso);
	putc(' ', stream);
	fprintf(stream, "%12.7g",  outobj->fluxerr_iso);
	putc(' ', stream);
	fprintf(stream, "%8.4f",  outobj->mag_iso);
	putc(' ', stream);
	fprintf(stream, "%8.4f",  outobj->magerr_iso);
	putc(' ', stream);

	fprintf(stream, "%12.7g",  outobj->flux_isocor);
	putc(' ', stream);
	fprintf(stream, "%12.7g",  outobj->fluxerr_isocor);
	putc(' ', stream);
	fprintf(stream, "%12.7g",  outobj->flux_aper[0]);
	putc(' ', stream);
	fprintf(stream, "%12.7g",  outobj->fluxerr_aper[0]);
	putc(' ', stream);
	fprintf(stream, "%8.4f",  outobj->mag_aper[0]);
	putc(' ', stream);
	fprintf(stream, "%8.4f",  outobj->magerr_aper[0]);
	putc(' ', stream);

	fprintf(stream, "%12.7g",  outobj->flux_auto);
	putc(' ', stream);
	fprintf(stream, "%12.7g",  outobj->fluxerr_auto);
	putc(' ', stream);
	fprintf(stream, "%8.4f",  outobj->mag_auto);
	putc(' ', stream);
	fprintf(stream, "%8.4f",  outobj->magerr_auto);
	putc(' ', stream);
	fprintf(stream, "%5.1f",  outobj->theta);
	putc(' ', stream);

	fprintf(stream, "%15.10e",  outobj->poserr_mx2);
	putc(' ', stream);
	fprintf(stream, "%15.10e",  outobj->poserr_my2);
	putc(' ', stream);
	fprintf(stream, "%15.10e",  outobj->poserr_mxy);
	putc(' ', stream);

	fprintf(stream, "%3d",  	outobj->flag);
	putc(' ', stream);
	fprintf(stream, "%8.3f",  	outobj->elong);
	putc(' ', stream);
	fprintf(stream, "%8.3f",  	outobj->ellip);
	putc(' ', stream);
	fprintf(stream, "%5.2f",  	outobj->sprob);
	putc(' ', stream);

	fprintf(stream, "%12.7g",  	outobj->peak);
	putc(' ', stream);
	fprintf(stream, "%5.2f",  	outobj->kronfactor);
	putc(' ', stream);
	fprintf(stream, "%12.7g",  	outobj->bkg);
	putc(' ', stream);
	fprintf(stream, "%8.2f",  	outobj->fwhm);
	//putc(' ', stream);

	putc('\n', stream);


	return;
}

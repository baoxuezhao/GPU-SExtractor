/*
 * cudaanalyse.h
 *
 *  Created on: 17 Mar, 2013
 *      Author: zhao
 */

#ifndef CUDAANALYSE_H_
#define CUDAANALYSE_H_

#include "../types.h"

void analyse_gpu(
		double		flagobj_poserr_mx2,
		int			flagobj_peakx,
		int			flagobj_iso_0,
		float		flagobj_fwhm,
		/* field global vars */
		float		backsig,
		double		gain,
		double		ngamma,
		float		thresh,
		float		dthresh,
		double		satur_level,
		/*	constant vars */
		int			plistexist_var,
		int			plistoff_var,
		int			prefs_weightgain_flag,
		int			prefs_ext_minarea,
		int			prefs_clean_flag,
		int			detect_type,
		int 		nbackx,
		int 		nbacky,
		int			backw,
		int			backh,
		int			width,
		int 		numobj);

void endobject_gpu(
		prefstruct	*h_prefs,
		brainstruct *h_brain,
		obj2struct	*h_flagobj2,
		outobjstruct *h_outobj,
		wcsstruct	*h_wcs,
		/*---field vars---*/
		int			fymin,
		int			fymax,
		int			width,
		int			height,
		float		backsig,
		float		thresh,
		double		satur_level,
		double		pixscale,
		double		ngamma,
		double		gain,
		int numobj);

//void copyobjects(objstruct *h_objectlist, int total);
void copyobjects(
		unsigned int *h_index,
		unsigned int *h_flag,
		unsigned int *h_fdnpix,
		float 		 *h_dthresh,

		unsigned int *h_xmin,
		unsigned int *h_xmax,
		unsigned int *h_ymin ,
		unsigned int *h_ymax,
		unsigned int *h_dnpix,
		unsigned int *h_npix,
		unsigned int *h_peakx,
		unsigned int *h_peaky ,

		float 	*h_bkg,
		float 	*h_dbkg ,
		float	*h_sigbkg,
		float 	*h_a,
		float 	*h_b,
		float 	*h_cxx,
		float 	*h_cxy,
		float 	*h_cyy,
		float 	*h_theta,
		float 	*h_abcor,
		float	*h_peak,
		float 	*h_dpeak,
		float 	*h_fdpeak ,
		float	*h_flux ,
		float 	*h_dflux,
		float 	*h_fdflux,
		float	*h_fluxerr,
		float	*h_thresh,
		float	*h_mthresh,
		float	*h_fwhm,

		double 	*h_mx,
		double 	*h_my,
		double 	*h_mx2,
		double 	*h_my2,
		double 	*h_mxy,
		double	*h_poserr_mx2 ,
		double	*h_poserr_my2,
		double	*h_poserr_mxy,
		char  	*h_singuflag,
		int		**h_iso,
		int total);

#endif /* CUDAANALYSE_H_ */

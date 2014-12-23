/*
 * cudaphotom.h
 *
 *  Created on: Apr 4, 2013
 *      Author: zhao
 */

#ifndef CUDAPHOTOM_H_
#define CUDAPHOTOM_H_

#include	"../types.h"

__device__ void  computeisocorflux1(
		unsigned int 	*d_fdnpix,
		float			*d_dthresh,
		float 			*d_flux,
		float			*d_fluxerr,
		float			*d_flux_isocor,
		float			*d_fluxerr_isocor,
		float			flagobj2_fluxerr_isocor,
		int tid);

__device__ void  computeaperflux1(

		float		*d_pixelArray,
		int			*d_labelArray, //allocation map
		unsigned int *d_flag,
		float		*d_dbkg,
		float		**d_flux_aper,
		float		**d_fluxerr_aper,
		float		**d_mag_aper,
		float		**d_magerr_aper,
		double		*d_mx,
		double		*d_my,
		prefstruct	*d_prefs,
		int 		width,
		int			height,
		int			fymin,
		int			fymax,
		double		ngamma,
		double		gain,
		float		backsig,
		int			tid,
		int 		i);

#endif /* CUDAPHOTOM_H_ */

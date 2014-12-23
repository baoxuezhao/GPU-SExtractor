/*
 				fitswcs.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	LDACTools+
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for fitswcs.c
*
*	Last modify:	26/04/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSWCS_H_
#define _FITSWCS_H_

#include "types.h"
/*-------------------------------- macros -----------------------------------*/

/*-------------------------------- typedefs ---------------------------------*/

/*------------------------------- structures --------------------------------*/

/*------------------------------- functions ---------------------------------*/

extern wcsstruct	*create_wcs(char **ctype, double *crval, double *crpix,
				double *cdelt, int *naxisn, int naxis),
			*copy_wcs(wcsstruct *wcsin),
			*read_wcs(tabstruct *tab);

extern double		fmod_0_p360(double angle),
			fmod_m90_p90(double angle),
			sextodegal(char *hms),
			sextodegde(char *dms),
			wcs_dist(wcsstruct *wcs,
				double *wcspos1, double *wcspos2),
			wcs_jacobian(wcsstruct *wcs, double *pixpos,
				double *jacob),
			wcs_scale(wcsstruct *wcs, double *pixpos);

extern int		celsys_to_eq(wcsstruct *wcs, double *wcspos),
			eq_to_celsys(wcsstruct *wcs, double *wcspos),
			raw_to_red(wcsstruct *wcs,
				double *pixpos, double *redpos),
			raw_to_wcs(wcsstruct *wcs,
				double *pixpos, double *wcspos),
			reaxe_wcs(wcsstruct *wcs, int lng, int lat),
			red_to_raw(wcsstruct *wcs,
				double *redpos, double *pixpos),
			wcs_chirality(wcsstruct *wcs),
			wcs_supproj(char *name),
			wcs_to_raw(wcsstruct *wcs,
				double *wcspos, double *pixpos);

extern char		*degtosexal(double alpha, char *str),
			*degtosexde(double delta, char *str);

extern void		b2j(double yearobs, double alphain, double deltain,
				double *alphaout, double *deltaout),
			end_wcs(wcsstruct *wcs),
			init_wcs(wcsstruct *wcs),
			init_wcscelsys(wcsstruct *wcs),
			invert_wcs(wcsstruct *wcs),
			frame_wcs(wcsstruct *wcsin, wcsstruct *wcsout),
			j2b(double yearobs, double alphain, double deltain,
				double *alphaout, double *deltaout),
			precess(double yearin, double alphain, double deltain,
				double yearout,
				double *alphaout, double *deltaout),
			precess_wcs(wcsstruct *wcs, double yearin,
				double yearout),
			range_wcs(wcsstruct *wcs),
			write_wcs(tabstruct *tab, wcsstruct *wcs);

#endif

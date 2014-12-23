/*
 				tnx.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	WCSlib
*
*	Author:		E.BERTIN (IAP), based on D.Mink (SAO) WCSTools
*
*	Contents:       Include to handle TNX astrometric format (from IRAF).
*
*
*	Last modify:	28/11/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _TNX_H_
#define _TNX_H_
#include "../types.h"

/*-------------------------------- macros -----------------------------------*/

#define		TNX_MAXCHARS	2048	/* Maximum FITS "WAT" string length */

/* TNX permitted types of surfaces */
#define		TNX_CHEBYSHEV	1
#define		TNX_LEGENDRE	2
#define		TNX_POLYNOMIAL	3

/* TNX cross-terms flags */
#define		TNX_XNONE	0	/* no x-terms (old no) */
#define		TNX_XFULL	1	/* full x-terms (new yes) */
#define		TNX_XHALF	2	/* half x-terms (new) */

/*----------------------------- Internal constants --------------------------*/

/*------------------------------- structures --------------------------------*/

/*------------------------------- functions ---------------------------------*/

tnxaxisstruct	*copy_tnxaxis(tnxaxisstruct *axis),
		*read_tnxaxis(char *tnxstr);

double		raw_to_tnxaxis(tnxaxisstruct *axis, double x, double y);

void		free_tnxaxis(tnxaxisstruct *axis);

#endif


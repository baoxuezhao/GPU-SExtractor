 /*
 				pattern.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Authors:	E.BERTIN (IAP)
*
*	Contents:	Include file for pattern.c.
*
*	Last modify:	19/11/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PROFIT_H_
#include "profit.h"
#endif

#ifndef _PATTERN_H_
#define _PATTERN_H_

#include 	"types.h"

/*-------------------------------- flags ------------------------------------*/
/*-------------------------------- macros -----------------------------------*/
/*----------------------------- Internal constants --------------------------*/

#define	PATTERN_FMAX	4	/* Maximum pattern angular frequency */
#define	PATTERN_NCOMP	16	/* Default number of components (radii) */
#define	PATTERN_SCALE	5.0	/* Pattern scale in units of r_eff */
#define	PATTERN_MARGIN	0.2	/* Pattern margin in fractions of radius */
#define	PATTERN_BTMAX	0.6	/* Maximum B/T for pure disk scaling */

/* NOTES:
One must have:	PATTERN_SIZE > 1
		PATTERN_SCALE > 0.0
		PATTERN_BTMAX < 1.0
*/

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

patternstruct	*pattern_init(profitstruct *profit, pattypenum ptype, int nvec);

float		pattern_spiral(patternstruct *pattern);
void		pattern_compmodarg(patternstruct *pattern,profitstruct *profit),
		pattern_create(patternstruct *pattern, profitstruct *profit),
		pattern_end(patternstruct *pattern),
		pattern_fit(patternstruct *pattern, profitstruct *profit);

#endif


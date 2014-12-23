 /*
 				photom.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP & Leiden observatory
*
*	Contents:	Include file for photom.h.
*
*	Last modify:	22/10/2004
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

/* NOTES:
One must have:	APER_OVERSAMP >= 1
		KRON_NSIG > 0.0
		PETRO_NSIG > 0.0
		CROWD_THRESHOLD >= 0
*/

/*------------------------------- functions ---------------------------------*/
extern void	computeaperflux(picstruct *, picstruct *, objstruct *, int),
		computeautoflux(picstruct *, picstruct *, picstruct *,
			picstruct *, objstruct *),
		computeisocorflux(picstruct *, objstruct *),
		computemags(picstruct *, objstruct *),
		computepetroflux(picstruct *, picstruct *, picstruct *,
				picstruct *, objstruct *);

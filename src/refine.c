/*
 refine.c

 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	Part of:	SExtractor
 *
 *	Author:		E.BERTIN (IAP)
 *
 *	Contents:	functions to refine extraction of objects.
 *
 *	Last modify:	10/01/2008
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"plist.h"
#include	"extract.h"

#ifndef	RAND_MAX
#define	RAND_MAX	2147483647
#endif
#define	NSONMAX			1024	/* max. number per level */
#define	NBRANCH			16	/* starting number per branch */


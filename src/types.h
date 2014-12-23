 /*
 				types.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global type definitions.
*
*	Last modify:	18/03/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
#ifndef __TYPES_H__
#define __TYPES_H__

#include <stdio.h>
#include "define.h"

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif
//#ifndef _FITSWCS_H_
//#include "fitswcs.h"
//#endif

/*--------------------------------- typedefs --------------------------------*/
typedef	unsigned char	BYTE;			/* a byte */
typedef	unsigned short	USHORT;			/* 0 to 65535 integers */
typedef	char			pliststruct;	/* Dummy type for plist */

typedef	int				LONG;
typedef	unsigned int	ULONG;

typedef  enum {BACK_RELATIVE, BACK_ABSOLUTE}
		backenum;				/* BACK_TYPE */

typedef  enum {CHECK_NONE, CHECK_IDENTICAL, CHECK_BACKGROUND,
        CHECK_BACKRMS, CHECK_MINIBACKGROUND, CHECK_MINIBACKRMS,
        CHECK_SUBTRACTED, CHECK_FILTERED, CHECK_OBJECTS, CHECK_APERTURES,
	CHECK_SEGMENTATION, CHECK_ASSOC, CHECK_SUBOBJECTS,
	CHECK_SUBPSFPROTOS, CHECK_PSFPROTOS,
	CHECK_SUBPCPROTOS, CHECK_PCPROTOS, CHECK_PCOPROTOS,
	CHECK_MAPSOM, CHECK_SUBPROFILES, CHECK_PROFILES, CHECK_PATTERNS,
	MAXCHECK}
		checkenum;				/* CHECK_IMAGE type */

typedef  enum {WEIGHT_NONE, WEIGHT_FROMBACK, WEIGHT_FROMRMSMAP,
		WEIGHT_FROMVARMAP, WEIGHT_FROMWEIGHTMAP, WEIGHT_FROMINTERP}
		weightenum;				/* WEIGHT_IMAGE type */

typedef  enum {CELSYS_NATIVE, CELSYS_PIXEL, CELSYS_EQUATORIAL, CELSYS_GALACTIC,
		  	CELSYS_ECLIPTIC, CELSYS_SUPERGALACTIC}
		celsysenum;
/*--------------------------------- objects ---------------------------------*/
/* I: "PIXEL" parameters */

typedef struct	__attribute__((aligned(16)))
  {
/* ---- basic parameters */
  int		number;				/* ID */
  int		fdnpix;				/* nb of extracted pix */
  int		dnpix;				/* nb of pix above thresh  */
  int		npix;				/* "" in measured frame */
  float		fdflux;				/* integrated ext. flux */
  float		dflux;				/* integrated det. flux */
  float		flux;				/* integrated mes. flux */
  float		fluxerr;			/* integrated variance */
  PIXTYPE	fdpeak;				/* peak intensity (ADU) */
  PIXTYPE	dpeak;				/* peak intensity (ADU) */
  PIXTYPE	peak;				/* peak intensity (ADU) */
/* ---- astrometric data */
  int		peakx,peaky;			/* pos of brightest pix */
  double       	mx, my;				/* barycenter */
  double	poserr_mx2, poserr_my2,
		poserr_mxy;			/* Error ellips moments */
/* ---- morphological data */			
  int		xmin,xmax,ymin,ymax,ycmin,ycmax;/* x,y limits */
  PIXTYPE	*blank, *dblank; 	       	/* BLANKing sub-images  */
  int		*submap;			/* Pixel-index sub-map */
  int		subx,suby, subw,subh;		/* sub-image pos. and size */
  short		flag;				/* extraction flags */
  BYTE		wflag;				/* weighted extraction flags */
  FLAGTYPE	imaflag[MAXFLAG];		/* flags from FLAG-images */
  BYTE		singuflag;			/* flags for singularities */
  int		imanflag[MAXFLAG];     		/* number of MOST flags */
  double	mx2,my2,mxy;			/* variances and covariance */
  float		a, b, theta, abcor;		/* moments and angle */
  float		cxx,cyy,cxy;			/* ellipse parameters */
  int		firstpix;			/* ptr to first pixel */
  int		lastpix;			/* ptr to last pixel */
  float		bkg, dbkg, sigbkg;		/* Background stats (ADU) */
  float		thresh;				/* measur. threshold (ADU) */
  float		dthresh;		       	/* detect. threshold (ADU) */
  float		mthresh;		       	/* max. threshold (ADU) */
  int		iso[NISO];			/* isophotal areas */
  float		fwhm;				/* IMAGE FWHM */
  }	objstruct;

/* II: "BLIND" parameters */
typedef struct	__attribute__((aligned(16)))
  {
/* ---- photometric data */
  float		flux_iso;			/* ISO integrated flux */
  float		fluxerr_iso;		/* RMS error on ISO flux */
  float		mag_iso;			/* ISO mag */
  float		magerr_iso;			/* ISO mag uncertainty */
  float		flux_isocor;		/* ISOCOR integrated flux */
  float		fluxerr_isocor;		/* RMS error on ISOCOR flux */
  float		mag_isocor;			/* ISOCOR mag */
  float		magerr_isocor;		/* ISOCOR mag uncertainty */
  float		kronfactor;			/* kron parameter */
  float		flux_auto;			/* AUTO integrated flux */
  float		fluxerr_auto;		/* RMS error on AUTO flux */
  float		mag_auto;			/* AUTO mag */
  float		magerr_auto;		/* AUTO mag uncertainty */
  float		petrofactor;		/* kron parameter */
  float		flux_petro;			/* AUTO integrated flux */
  float		fluxerr_petro;		/* RMS error on AUTO flux */
  float		mag_petro;			/* AUTO mag */
  float		magerr_petro;		/* AUTO mag uncertainty */
  float		flux_best;			/* BEST integrated flux */
  float		fluxerr_best;		/* RMS error on BEST flux */
  float		mag_best;			/* BEST mag */
  float		magerr_best;		/* BEST mag uncertainty */
  float		flux_aper[NAPER];		/* APER flux vector */
  float		fluxerr_aper[NAPER];	/* APER flux error vector  */
  float		mag_aper[NAPER];		/* APER magnitude vector */
  float		magerr_aper[NAPER];		/* APER mag error vector */
  float		flux_win;			/* WINdowed flux*/
  float		fluxerr_win;		/* WINdowed flux error */
  float		mag_win;			/* WINdowed magnitude */
  float		magerr_win;			/* WINdowed magnitude error */
/* ---- astrometric data */
  double	posx,posy;				/* "FITS" pos. in pixels */
  double	jacob[NAXIS*NAXIS];		/* Local deproject. Jacobian */
  double	mamaposx,mamaposy;		/* "MAMA" pos. in pixels */
  float		sposx,sposy;			/* single precision pos. */
  float		poserr_a, poserr_b,
		poserr_theta;			/* Error ellips parameters */
  float		poserr_cxx, poserr_cyy,
		poserr_cxy;			/* pos. error ellipse */
  double	poserr_mx2w, poserr_my2w,
		poserr_mxyw;			/* WORLD error moments */
  float		poserr_aw, poserr_bw,
		poserr_thetaw;			/* WORLD error parameters */
  float		poserr_thetas;			/* native error pos. angle */
  float		poserr_theta2000;		/* J2000 error pos. angle */
  float		poserr_theta1950;		/* B1950 error pos. angle */
  float		poserr_cxxw, poserr_cyyw,
		poserr_cxyw;			/* WORLD error ellipse */
  double	mx2w,my2w,mxyw;			/* WORLD var. and covar. */
  double	peakxw, peakyw;			/* WORLD of brightest pix */
  double	mxw, myw;			/* WORLD barycenters */
  double	alphas, deltas;			/* native alpha, delta */
  float		thetas;				/* native position angle E/N*/
  double	peakalphas, peakdeltas;		/* native for brightest pix */
  double	peakalpha2000, peakdelta2000;	/* J2000 for brightest pix */
  double	peakalpha1950, peakdelta1950;	/* B1950 for brightest pix */
  double	alpha2000, delta2000;		/* J2000 alpha, delta */
  float		theta2000;			/* J2000 position angle E/N */
  double	dtheta2000;			/* North J2000 - native angle*/
  double	alpha1950, delta1950;		/* B1950 alpha, delta */
  float		theta1950;			/* B1950 position angle E/N */
  double	dtheta1950;			/* North B1950 - native angle*/
  float		aw, bw;				/* WORLD ellipse size */
  float		thetaw;				/* WORLD position angle */
  float		cxxw,cyyw,cxyw;			/* WORLD ellipse parameters */
  float		npixw, fdnpixw;			/* WORLD isophotal areas */
  float		threshmu;			/* det. surface brightnees */
  float		maxmu;				/* max. surface brightnees */
  float		elong;				/* elongation */
  float		ellip;				/* ellipticity */
  float		polar;				/* Kaiser's "polarization" */
  float		polarw;				/* WORLD "polarization" */
  float		sprob;				/* Stellarity index */
  float		fwhmw;				/* WORLD FWHM */
  float		*assoc;				/* ASSOCiated data */
  int		assoc_number;			/* nb of ASSOCiated objects */
  float		*vignet;			/* Pixel data */
  float		*vigshift;			/* (Shifted) pixel data */

///* Windowed measurements */
//  double	winpos_x,winpos_y;		/* Windowed barycenter */
//  double	winposerr_mx2, winposerr_my2,
//		winposerr_mxy;			/* Error ellips moments */
//  float		winposerr_a, winposerr_b,
//		winposerr_theta;		/* Error ellips parameters */
//  float		winposerr_cxx, winposerr_cyy,
//		winposerr_cxy;			/* pos. error ellipse */
//  double	winposerr_mx2w, winposerr_my2w,
//		winposerr_mxyw;			/* WORLD error moments */
//  float		winposerr_aw, winposerr_bw,
//		winposerr_thetaw;		/* WORLD error parameters */
//  float		winposerr_thetas;		/* native error pos. angle */
//  float		winposerr_theta2000;		/* J2000 error pos. angle */
//  float		winposerr_theta1950;		/* B1950 error pos. angle */
//  float		winposerr_cxxw, winposerr_cyyw,
//		winposerr_cxyw;			/* WORLD error ellipse */
//  double	win_mx2, win_my2,
//		win_mxy;			/* Windowed moments */
//  float		win_a, win_b,
//		win_theta;			/* Windowed ellipse parameters*/
//  float		win_polar;			/* Windowed "polarization" */
//  float		win_cxx, win_cyy,
//		win_cxy;			/* Windowed ellipse parameters*/
//  double	win_mx2w, win_my2w,
//		win_mxyw;			/* WORLD windowed moments */
//  float		win_aw, win_bw,
//		win_thetaw;			/* WORLD ellipse parameters */
//  float		win_polarw;			/* WORLD WIN "polarization" */
//  float		win_thetas;		/* native error pos. angle */
//  float		win_theta2000;		/* J2000 error pos. angle */
//  float		win_theta1950;		/* B1950 error pos. angle */
//  float		win_cxxw, win_cyyw,
//		win_cxyw;			/* WORLD ellipse parameters */
//  double	winpos_xw, winpos_yw;		/* WORLD coordinates */
//  double	winpos_alphas, winpos_deltas;	/* native alpha, delta */
//  double	winpos_alpha2000, winpos_delta2000;	/* J2000 alpha, delta */
//  double	winpos_alpha1950, winpos_delta1950;	/* B1950 alpha, delta */
//  short		winpos_niter;			/* Number of WIN iterations */
//  short		win_flag;			/* 1:x2<0 2:xy=x2 4:flux<0 */
//
// /* ---- SOM fitting */
//  float		flux_somfit;			/* Fitted amplitude */
//  float		fluxerr_somfit;			/* RMS error on SOM flux */
//  float		mag_somfit;			/* Magnitude from SOM fit */
//  float		magerr_somfit;			/* Mag. err. from SOM fit */
//  float		stderr_somfit;			/* Fitting reduced error */
//  float		*vector_somfit;			/* SOM fit vector position */
///* ---- Growth curves and stuff */
//  float		*flux_growth;			/* Cumulated growth_curve */
//  float		flux_growthstep;		/* Growth-curve step */
//  float		*mag_growth;			/* Cumulated growth_curve */
//  float		mag_growthstep;			/* Growth-curve step */
//  float		*flux_radius;			/* f-light-radii */
//  float		hl_radius;			/* Scalar half-light radius */
///* ---- PSF-fitting */
//  float		*flux_psf;			/* Flux from PSF-fitting */
//  float		*fluxerr_psf;			/* RMS error on PSF flux */
//  float		*mag_psf;			/* Mag from PSF-fitting */
//  float		*magerr_psf;			/* RMS mag from PSF-fitting */
//  float		*x_psf, *y_psf;			/* Coords from PSF-fitting */
//  short		niter_psf;			/* # of PSF-fitting iterat. */
//  short		npsf;				/* # of fitted PSFs */
//  float		chi2_psf;			/* Red. chi2 of PSF-fitting */
//  double	xw_psf, yw_psf;			/* WORLD coords */
//  double	alphas_psf, deltas_psf;		/* native alpha, delta */
//  double	alpha2000_psf, delta2000_psf;	/* J2000 alpha, delta */
//  double	alpha1950_psf, delta1950_psf;	/* B1950 alpha, delta */
//  double	poserrmx2_psf, poserrmy2_psf,
//		poserrmxy_psf;			/* Error ellips moments */
//  float		poserra_psf, poserrb_psf,
//		poserrtheta_psf;		/* Error ellips parameters */
//  float		poserrcxx_psf, poserrcyy_psf,
//		poserrcxy_psf;			/* pos. error ellipse */
//  double	poserrmx2w_psf, poserrmy2w_psf,
//		poserrmxyw_psf;			/* WORLD error moments */
//  float		poserraw_psf, poserrbw_psf,
//		poserrthetaw_psf;		/* WORLD error parameters */
//  float		poserrthetas_psf;		/* native error pos. angle */
//  float		poserrtheta2000_psf;		/* J2000 error pos. angle */
//  float		poserrtheta1950_psf;		/* B1950 error pos. angle */
//  float		poserrcxxw_psf, poserrcyyw_psf,
//		poserrcxyw_psf;			/* WORLD error ellipse */
///* ---- PC-fitting */
//  double	mx2_pc,my2_pc,mxy_pc;		/* PC 2nd-order parameters */
//  float		a_pc,b_pc,theta_pc;		/* PC shape parameters */
//  float		*vector_pc;			/* Principal components */
//  float		gdposang;			/* Gal. disk position angle */
//  float		gdscale;			/* Gal. disk scalelength */
//  float		gdaspect;			/* Gal. disk aspect-ratio */
//  float		gde1,gde2;			/* Gal. disk ellipticities */
//  float		gbratio;			/* Galaxy B/T */
//  float		gbposang;			/* Gal. bulge position angle */
//  float		gbscale;			/* Gal. bulge scalelength */
//  float		gbaspect;			/* Gal. bulge aspect-ratio */
//  float		gbe1,gbe2;			/* Gal. bulge ellipticities */
//  float		flux_galfit;			/* Galaxy tot. flux from fit */
//  float		fluxerr_galfit;			/* RMS error on galfit flux */
//  float		mag_galfit;			/* Galaxy tot. mag from fit */
//  float		magerr_galfit;			/* RMS error on galfit mag */
///* ---- Profile-fitting */
//  float		*prof_vector;			/* Profile parameters */
//  float		*prof_errvector;		/* Profile parameter errors */
//  float		prof_chi2;			/* Reduced chi2 */
//  BYTE		prof_flag;			/* Model-fitting flags */
//  BYTE		prof_flagw;			/* Model-fitting WORLD flag */
//  short		prof_niter;			/* # of model-fitting iter. */
//  float		flux_prof;			/* Flux from model-fitting */
//  float		fluxerr_prof;			/* RMS error on model flux */
//  float		mag_prof;			/* Mag from model-fitting */
//  float		magerr_prof;			/* RMS mag from model-fitting */
//  float		x_prof, y_prof;			/* Coords from model-fitting*/
//  double	xw_prof, yw_prof;		/* WORLD coords */
//  double	alphas_prof, deltas_prof;	/* native alpha, delta */
//  double	alpha2000_prof, delta2000_prof;	/* J2000 alpha, delta */
//  double	alpha1950_prof, delta1950_prof;	/* B1950 alpha, delta */
//  double	poserrmx2_prof, poserrmy2_prof,
//		poserrmxy_prof;			/* Error ellips moments */
//  float		poserra_prof, poserrb_prof,
//		poserrtheta_prof;		/* Error ellips parameters */
//  float		poserrcxx_prof, poserrcyy_prof,
//		poserrcxy_prof;			/* pos. error ellipse */
//  double	poserrmx2w_prof, poserrmy2w_prof,
//		poserrmxyw_prof;		/* WORLD error moments */
//  float		poserraw_prof, poserrbw_prof,
//		poserrthetaw_prof;		/* WORLD error parameters */
//  float		poserrthetas_prof;		/* native error pos. angle */
//  float		poserrtheta2000_prof;		/* J2000 error pos. angle */
//  float		poserrtheta1950_prof;		/* B1950 error pos. angle */
//  float		poserrcxxw_prof, poserrcyyw_prof,
//		poserrcxyw_prof;		/* WORLD error ellipse */
//  double	prof_mx2, prof_my2, prof_mxy;	/* Profile model moments */
//  double	prof_mx2w, prof_my2w, prof_mxyw;/* WORLD profile model moments*/
//  float		prof_eps1, prof_eps2;		/* Profile model ellip.vector */
//  float		prof_e1, prof_e2;		/* Profile model ellip.vector */
//  float		prof_class_star;		/* Model-fitting star/gal class*/
//  float		prof_concentration;		/* Model-fitting concentration*/
//  float		prof_offset_flux;		/* Background offset */
//  float		prof_offset_fluxerr;		/* RMS error */
//  float		prof_spheroid_flux;		/* Spheroid total flux */
//  float		prof_spheroid_fluxerr;		/* RMS error */
//  float		prof_spheroid_mag;		/* Spheroid "total" mag */
//  float		prof_spheroid_magerr;		/* RMS error */
//  float		prof_spheroid_reff;		/* Spheroid effective radius */
//  float		prof_spheroid_refferr;		/* RMS error */
//  float		prof_spheroid_reffw;		/* WORLD spheroid eff. radius */
//  float		prof_spheroid_refferrw;		/* RMS error */
//  float		prof_spheroid_aspect;		/* Spheroid aspect ratio */
//  float		prof_spheroid_aspecterr;	/* RMS error */
//  float		prof_spheroid_aspectw;		/* WORLD spheroid aspect ratio*/
//  float		prof_spheroid_aspecterrw;	/* RMS error */
//  float		prof_spheroid_theta;		/* Spheroid position angle */
//  float		prof_spheroid_thetaerr;		/* RMS error */
//  float		prof_spheroid_thetaw;		/* WORLD spheroid pos. angle */
//  float		prof_spheroid_thetaerrw;	/* RMS error */
//  float		prof_spheroid_thetas;		/* Sky spheroid pos. angle */
//  float		prof_spheroid_theta2000;	/* J2000 spheroid pos. angle */
//  float		prof_spheroid_theta1950;	/* B1950 spheroid pos. angle */
//  float		prof_spheroid_sersicn;		/* Spheroid Sersic index */
//  float		prof_spheroid_sersicnerr;	/* RMS error */
//  float		prof_disk_flux;			/* Disk total flux */
//  float		prof_disk_fluxerr;		/* RMS error */
//  float		prof_disk_mag;			/* Disk "total" mag */
//  float		prof_disk_magerr;		/* RMS error */
//  float		prof_disk_scale;		/* Disk scale length */
//  float		prof_disk_scaleerr;		/* RMS error */
//  float		prof_disk_scalew;		/* WORLD disk scale length */
//  float		prof_disk_scaleerrw;		/* RMS error */
//  float		prof_disk_aspect;		/* Disk aspect ratio */
//  float		prof_disk_aspecterr;		/* RMS error */
//  float		prof_disk_aspectw;		/* WORLD disk aspect ratio */
//  float		prof_disk_aspecterrw;		/* RMS error */
//  float		prof_disk_inclination;		/* Disk inclination */
//  float		prof_disk_inclinationerr;	/* RMS error */
//  float		prof_disk_theta;		/* Disk position angle */
//  float		prof_disk_thetaerr;		/* RMS error */
//  float		prof_disk_thetaw;		/* WORLD disk position angle */
//  float		prof_disk_thetaerrw;		/* RMS error */
//  float		prof_disk_thetas;		/* Sky disk position angle */
//  float		prof_disk_theta2000;		/* J2000 disk position angle */
//  float		prof_disk_theta1950;		/* B1950 disk position angle */
//  float		*prof_disk_patternvector;	/* Disk pattern coefficients */
//  float		*prof_disk_patternmodvector;	/* Disk pattern moduli */
//  float		*prof_disk_patternargvector;	/* Disk pattern arguments */
//  float		prof_disk_patternspiral;	/* Disk pattern spiral index */
//  float		prof_bar_flux;			/* Galactic bar total flux */
//  float		prof_bar_fluxerr;		/* RMS error */
//  float		prof_bar_mag;			/* Bar "total" magnitude */
//  float		prof_bar_magerr;		/* RMS error */
//  float		prof_bar_length;		/* Bar length */
//  float		prof_bar_lengtherr;		/* RMS error */
//  float		prof_bar_lengthw;		/* WORLD bar length */
//  float		prof_bar_lengtherrw;		/* RMS error */
//  float		prof_bar_aspect;		/* Bar aspect ratio */
//  float		prof_bar_aspecterr;		/* RMS error */
//  float		prof_bar_aspectw;		/* WORLD bar aspect ratio */
//  float		prof_bar_aspecterrw;		/* RMS error */
//  float		prof_bar_posang;		/* Bar true prosition angle */
//  float		prof_bar_posangerr;		/* RMS error */
//  float		prof_bar_theta;			/* Bar projected angle */
//  float		prof_bar_thetaerr;		/* RMS error */
//  float		prof_bar_thetaw;		/* WORLD bar projected angle */
//  float		prof_bar_thetaerrw;		/* RMS error */
//  float		prof_bar_thetas;		/* Sky bar projected angle */
//  float		prof_bar_theta2000;		/* J2000 bar projected angle */
//  float		prof_bar_theta1950;		/* B1950 bar projected angle */
//  float		prof_arms_flux;			/* Spiral arms total flux */
//  float		prof_arms_fluxerr;		/* RMS error */
//  float		prof_arms_mag;			/* Arms "total" magnitude */
//  float		prof_arms_magerr;		/* RMS error */
//  float		prof_arms_scale;		/* Arms scalelength */
//  float		prof_arms_scaleerr;		/* RMS error */
//  float		prof_arms_scalew;		/* WORLD arms scalelength */
//  float		prof_arms_scaleerrw;		/* RMS error */
//  float		prof_arms_posang;		/* Arms true position angle */
//  float		prof_arms_posangerr;		/* RMS error */
//  float		prof_arms_thetaw;		/* WORLD arms position angle */
//  float		prof_arms_thetas;		/* Sky arms position angle */
//  float		prof_arms_theta2000;		/* J2000 arms position angle */
//  float		prof_arms_theta1950;		/* B1950 arms position angle */
//  float		prof_arms_pitch;		/* Arms pitch angle */
//  float		prof_arms_pitcherr;		/* RMS error */
//  float		prof_arms_start;		/* Arms starting radius */
//  float		prof_arms_starterr;		/* RMS error */
//  float		prof_arms_startw;		/* WORLD arms starting radius */
//  float		prof_arms_starterrw;		/* RMS error */
//  float		prof_arms_quadfrac;		/* Arms quadrature fraction */
//  float		prof_arms_quadfracerr;		/* RMS error */
/* ---- MEF */
  short		ext_number;			/* FITS extension number */
  }	obj2struct;

typedef struct __attribute__((aligned(16)))
{
	int 	number;
	double 	alpha2000;
	double 	delta2000;
	float 	sposx;
	float	sposy;
	float	flux_iso;
	float	fluxerr_iso;
	float	mag_iso;
	float	magerr_iso;
	float	flux_isocor;
	float	fluxerr_isocor;

	float	flux_aper[NAPER];
	float	fluxerr_aper[NAPER];
	float	mag_aper[NAPER];
	float	magerr_aper[NAPER];

	float	flux_auto;
	float	fluxerr_auto;
	float	mag_auto;
	float	magerr_auto;
	float	theta;

	double 	poserr_mx2;	//ERRX2_IMAGE
	double 	poserr_my2;	//ERRY2_IMAGE
	double	poserr_mxy;	//ERRXY_IMAGE

	short	flag;		//FLAGS
	float	elong;		//ELONGATION
	float	ellip;		//ELLIPTICITY
	float	sprob;		//CLASS_STAR
	float	peak;		//FLUX_MAX
	float	kronfactor;	//KRON_RADIUS
	float	bkg;		//BACKGROUND
	float	fwhm;		//FWHM_IMAGE

} outobjstruct;
/*----------------------------- lists of objects ----------------------------*/
typedef struct
  {
  int		nobj;			/* number of objects in list */
  objstruct	*obj;			/* pointer to the object array */
  int		npix;			/* number of pixels in pixel-list */
  pliststruct	*plist;			/* pointer to the pixel-list */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
  }	objliststruct;


typedef struct  __attribute__((aligned(16)))
{
     int flag;
     char pcode[4];
     char lngtyp[5], lattyp[5];
     int lng, lat;
     int cubeface;
} wcsprm;

typedef struct __attribute__((aligned(16)))
{
   int flag;
   int naxis;
   double *crpix;
   double *pc;
   double *cdelt;

   /* Intermediates. */
   double *piximg;
   double *imgpix;
} linprm;

typedef struct /*tnxaxis*/ __attribute__((aligned(16)))
{
  int		type;			/* Projection correction type */
  int		xorder,yorder;	/* Polynomial orders */
  int		xterms;			/* Well... */
  int		ncoeff;			/* Number of polynom coefficients */
  double	xrange,yrange;	/* Coordinate ranges */
  double	xmaxmin,ymaxmin;/* Well... */
  double	*coeff;			/* Polynom coefficients */
  double	*xbasis,*ybasis;/* Basis function values */
} tnxaxisstruct;

typedef struct /*poly*/ __attribute__((aligned(16)))
{
  double	*basis;		/* Current values of the basis functions */
  double	*coeff;		/* Polynom coefficients */
  int		ncoeff;		/* Number of coefficients */
  int		*group;		/* Groups */
  int		ndim;		/* dimensionality of the polynom */
  int		*degree;	/* Degree in each group */
  int		ngroup;		/* Number of different groups */
} polystruct;

typedef struct  __attribute__((aligned(16)))
{
   int flag;
   int n;
   double r0;
   double p[200];
   double w[10];
   tnxaxisstruct	*tnx_latcor;
   tnxaxisstruct	*tnx_lngcor;
   polystruct		*inv_x;
   polystruct		*inv_y;
} prjprm;

typedef struct  __attribute__((aligned(16)))
{
   int flag;
   double ref[4];
   double euler[5];

#if __STDC__  || defined(__cplusplus)
   int (*prjfwd)(const double,
		   	   	   const double,
		   	   	   prjprm *,
		   	   	   double *, double *);
   int (*prjrev)(const double,
		   	   	   const double,
		   	   	   prjprm *,
                 double *, double *);
#else
   int (*prjfwd)();
   int (*prjrev)();
#endif
} celprm;

typedef struct /*wcs*/ __attribute__((aligned(16)))
{
    int		naxis;				/* Number of image axes */
    int		naxisn[NAXIS];		/* FITS NAXISx parameters */
    char	ctype[NAXIS][9];	/* FITS CTYPE strings */
    char	cunit[NAXIS][32];	/* FITS CUNIT strings */
    double	crval[NAXIS];		/* FITS CRVAL parameters */
    double	cdelt[NAXIS];		/* FITS CDELT parameters */
    double	crpix[NAXIS];		/* FITS CRPIX parameters */
    double	crder[NAXIS];		/* FITS CRDER parameters */
    double	csyer[NAXIS];		/* FITS CSYER parameters */
    double	cd[NAXIS*NAXIS];	/* FITS CD matrix */
    double	*projp;				/* FITS PV/PROJP mapping parameters */
    int		nprojp;				/* number of useful projp parameters */
    double	longpole,latpole;	/* FITS LONGPOLE and LATPOLE */
    double	wcsmin[NAXIS];		/* minimum values of WCS coords */
    double	wcsmax[NAXIS];		/* maximum values of WCS coords */
    double	wcsscale[NAXIS];	/* typical pixel scale at center */
    double	wcsscalepos[NAXIS];	/* WCS coordinates of scaling point */
    double	wcsmaxradius;		/* Maximum distance to wcsscalepos */
    int		outmin[NAXIS];		/* minimum output pixel coordinate */
    int		outmax[NAXIS];		/* maximum output pixel coordinate */
    int		lat,lng;		/* longitude and latitude axes # */
    double	r0;				/* projection "radius" */
    double	lindet;			/* Determinant of the local matrix */
    int		chirality;		/* Chirality of the CD matrix */
    double	pixscale;		/* (Local) pixel scale */
    double	ap2000,dp2000;	/* J2000 coordinates of pole */
    double	ap1950,dp1950;	/* B1950 coordinates of pole */
    double	obsdate;		/* Date of observations */
    double	equinox;		/* Equinox of observations */
    double	epoch;			/* Epoch of observations (deprec.) */
    enum {RDSYS_ICRS, RDSYS_FK5, RDSYS_FK4, RDSYS_FK4_NO_E, RDSYS_GAPPT}
  		radecsys;			/* FITS RADECSYS reference frame */
    celsysenum	celsys;		/* Celestial coordinate system */
    double	celsysmat[4];	/* Equ. <=> Cel. system parameters */
    int		celsysconvflag;	/* Equ. <=> Cel. conversion needed? */
    wcsprm	*wcsprm;		/* WCSLIB's wcsprm structure */
    linprm	*lin;			/* WCSLIB's linprm structure */
    celprm	*cel;			/* WCSLIB's celprm structure */
    prjprm 	*prj;			/* WCSLIB's prjprm structure */
    tnxaxisstruct *tnx_latcor;	/* IRAF's TNX latitude corrections */
    tnxaxisstruct *tnx_lngcor;	/* IRAF's TNX longitude corrections */
    polystruct	*inv_x;			/* Proj. correction polynom in x */
    polystruct	*inv_y;			/* Proj. correction polynom in y */
}	wcsstruct;

/*----------------------------- image parameters ----------------------------*/
typedef struct pic
  {
  char		filename[MAXCHAR];	/* pointer to the image filename */
  char		*rfilename;		/* pointer to the reduced image name */
  char		ident[MAXCHAR];		/* field identifier (read from FITS)*/
  char		rident[MAXCHAR];	/* field identifier (relative) */
  catstruct	*cat;			/* FITS structure */
  tabstruct	*tab;			/* FITS extension structure */
  FILE		*file;			/* pointer the image file structure */
/* ---- main image parameters */
  int		bitpix, bytepix;	/* nb of bits and bytes per pixel */
  int		bitsgn;			/* non-zero if signed integer data */
  int		width, height;		/* x,y size of the field */
  KINGSIZE_T	npix;			/* total number of pixels */
  double	bscale, bzero;		/* FITS scale and offset */
  double	ngamma;			/* normalized photo gamma */
  int		nlevels;		/* nb of quantification levels */
  float		pixmin, pixmax;		/* min and max values in frame */
  int		y;			/* y current position in field */
  int		ymin;			/* y limit (lowest accessible) */
  int		ymax;			/* y limit (highest accessible+1) */
  int		yblank;			/* y blanking limit (highest+1) */
  PIXTYPE	*strip;			/* pointer to the image buffer */
  FLAGTYPE	*fstrip;		/* pointer to the FLAG buffer */
  int		stripheight;		/* height  of a strip (in lines) */
  int		stripmargin;		/* number of lines in margin */
  int		stripstep;		/* number of lines at each read */
  int		stripy;			/* y position in buffer */
  int		stripylim;		/* y limit in buffer */
  int		stripysclim;		/* y scroll limit in buffer */
/* ---- basic astrometric parameters */
   double	pixscale;		/* pixel size in arcsec.pix-1 */
   double	epoch;			/* epoch of coordinates */
/* ---- basic photometric parameters */
   double	gain;			/* conversion factor in e-/ADU */
   double	satur_level;		/* saturation level in ADUs */
/* ---- background parameters */
  float		*back;			/* ptr to the background map in mem */
  float		*dback;			/* ptr to the background deriv. map */
  float		*sigma;			/* ptr to the sigma map */
  float		*dsigma;		/* Ptr to the sigma deriv. map */
  int		backw, backh;		/* x,y size of a bkgnd mesh */
  int		nbackp;			/* total nb of pixels per bkgnd mesh */
  int		nbackx, nbacky;		/* x,y number of bkgnd meshes */
  int		nback;			/* total number of bkgnd meshes */
  int		nbackfx, nbackfy;	/* x,y size of bkgnd filtering mask */
  float		backmean;		/* median bkgnd value in image */
  float		backsig;		/* median bkgnd rms in image */
  float		sigfac;			/* scaling RMS factor (for WEIGHTs) */
  PIXTYPE	*backline;		/* current interpolated bkgnd line */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
  backenum	back_type;		/* Background type */
/* ---- astrometric parameters */
  wcsstruct	*wcs;			/* astrometric data */
  struct structassoc	*assoc;		/* ptr to the assoc-list */
  int		flags;			/* flags defining the field type */
/* ---- image interpolation */
  int		interp_flag;		/* interpolation for this field? */
  PIXTYPE	*interp_backup;		/* backup line for interpolation */
  PIXTYPE	weight_thresh;		/* interpolation threshold */
  int		*interp_ytimeoutbuf;	/* interpolation timeout line buffer */
  int		interp_xtimeout;	/* interpolation timeout value in x */
  int		interp_ytimeout;	/* interpolation timeout value in y */
  struct pic	*reffield;	       	/* pointer to a reference field */
  OFF_T		mefpos;			/* Position in a MEF file */
  }	picstruct;

  /*----------------------------- Internal constants --------------------------*/

  typedef enum	{PNONE, FIXED, AUTO} aperttype;
  typedef enum	{CAT_NONE, ASCII, ASCII_HEAD, ASCII_SKYCAT, ASCII_VO,
    	FITS_LDAC, FITS_TPX, FITS_10}	 cattype;

  typedef enum	{THRESH_RELATIVE, THRESH_ABSOLUTE} threshtype;
  //typedef enum	{CCD, PHOTO} 	detecttype;

  typedef enum	{FLAG_OR, FLAG_AND, FLAG_MIN, FLAG_MAX, FLAG_MOST} flagtype;
  typedef enum	{IMAGE, AFILE}			back_origintype;
  typedef enum	{GLOBAL, LOCAL}			pbacktype;
  typedef enum	{QUIET, NORM, WARN, FULL} verbosetype;
  typedef enum	{ASSOC_FIRST, ASSOC_NEAREST, ASSOC_MEAN, ASSOC_MAGMEAN,
  	 ASSOC_SUM, ASSOC_MAGSUM, ASSOC_MIN, ASSOC_MAX} assoctype;

  typedef enum	{ASSOCSELEC_ALL, ASSOCSELEC_MATCHED, ASSOCSELEC_NOMATCHED}
  assocselectype;

  //typedef enum 	{MASK_NONE, MASK_BLANK, MASK_CORRECT} masktype;
  typedef enum	{INTERP_NONE, INTERP_VARONLY, INTERP_ALL} interptype;
  typedef enum	{PSFDISPLAY_SPLIT, PSFDISPLAY_VECTOR} psfdisplaytype;

  /* NOTES:
  One must have:	MAXLIST >= 1 (preferably >= 16!)
  */

  /*--------------------------------- typedefs --------------------------------*/

  typedef enum		{PATTERN_QUADRUPOLE, PATTERN_OCTOPOLE,
  			PATTERN_POLARFOURIER, PATTERN_POLARSHAPELETS,
  			PATTERN_NPATTERNS}
  				pattypenum; /* Pattern code */

  /*--------------------------- structure definitions -------------------------*/

  typedef struct
    {
    pattypenum	type;		/* Pattern code */
    int		ncomp;			/* Number of independent components */
    int		nmodes;			/* Number of modes per component */
    int		nfreq;			/* Number of waves per component */
    double	x[2];			/* Coordinate vector */
    double	rmax;			/* Largest radius in units of scale */
    double	*r;			/* Reduced radius */
    double	*norm;			/* Pattern vector norm */
    double	*coeff;			/* Fitted pattern coefficients */
    double	*mcoeff;		/* Modulus from pattern coefficients */
    double	*acoeff;		/* Argument from pattern coefficients */
    double	*modpix;		/* Pattern pixmaps */
    PIXTYPE	*lmodpix;		/* Low resolution pattern pixmaps */
    int		size[3];		/* Pixmap size for each axis */
    }	patternstruct;

  /*------------------------------- preferences -------------------------------*/
  typedef struct  __attribute__((aligned(16)))
    {
    char		**command_line;			/* Command line */
    int			ncommand_line;			/* nb of params */
    char		prefs_name[MAXCHAR];	/* prefs filename*/
    char		*(image_name[2]);		/* image filenames */
    int			nimage_name;			/* nb of params */
    char		cat_name[MAXCHAR];		/* catalog filename*/
  /*----- thresholding */
    double		dthresh[2];				/* detect. threshold */
    int			ndthresh;				/* (1 or 2 entries) */
    double		thresh[2];				/* analysis thresh. */
    int			nthresh;				/* (1 or 2 entries) */
    threshtype	thresh_type[2];			/* bkgnd type */

    int			nthresh_type;			/* nb of params */
  /*----- extraction */
    int			dimage_flag;			/* detect. image ? */
    int			ext_minarea;			/* min area in pix. */
    int			deb_maxarea;			/* max deblend. area */
    int			filter_flag;			/* smoothing on/off */
    char		filter_name[MAXCHAR];	/* mask filename */
    double		filter_thresh[2];			/* Filter thresholds */
    int			nfilter_thresh;				/* nb of params */
    int			deblend_nthresh;			/* threshold number */
    double		deblend_mincont;			/* minimum contrast */
    char		satur_key[8];				/* saturation keyword */
    double		satur_level;				/* saturation level */
    enum	{CCD, PHOTO}	detect_type;	/* detection type */
  /*----- Flagging */
    char		*(fimage_name[MAXFLAG]);	/* flagmap filenames */
    int			nfimage_name;				/* nb of params */
    flagtype	flag_type[MAXFLAG];	/* flag combination */
    int			imaflag_size;				/* requested iso nb1 */
    int			imanflag_size;				/* requested iso nb2 */
    int			nimaisoflag;				/* effective iso nb */
    int			nimaflag;				/* effective ima nb */
  /*----- cleaning */
    int			clean_flag;				/* allow cleaning ? */
    double		clean_param;				/* cleaning effic. */
  /*----- Weighting */
    char		*(wimage_name[2]);       	/* weight filenames */
    int			nwimage_name;				/* nb of params */
    weightenum	weight_type[2];				/* weighting scheme */
    int			nweight_type;				/* nb of params */
    int			weight_flag;				/* do we weight ? */
    int			dweight_flag;				/* detection weight? */
    int			weightgain_flag;			/* weight gain? */
  /*----- photometry */
    cattype		cat_type;					/* type of catalog */
    aperttype	apert_type;	/* type of aperture */
    double		apert[MAXNAPER];			/* apert size (pix) */
    int			naper;						/* effective apert. */
    int			flux_apersize, fluxerr_apersize;	/* requested apert. */
    int			mag_apersize, magerr_apersize;		/* requested apert. */
    double		autoparam[2];				/* Kron parameters */
    int			nautoparam;					/* nb of Kron params */
    double		petroparam[2];				/* Kron parameters */
    int		npetroparam;					/* nb of Kron params */
    double	autoaper[2];					/* minimum apertures */
    int		nautoaper;						/* nb of min. aperts */
    double	mag_zeropoint;					/* magnitude offsets */
    double	mag_gamma;						/* for emulsions */
    double	gain;							/* only for CCD */
    char		gain_key[8];				/* gain keyword */
  /*----- S/G separation */
    double	pixel_scale;					/* in arcsec */
    double	seeing_fwhm;					/* in arcsec */
    char		nnw_name[MAXCHAR];			/* nnw filename */
  /*----- background */
    back_origintype			back_origin;	/* origin of bkgnd */
    char		back_name[MAXCHAR];			/* bkgnd filename */
    backenum	back_type[2];				/* bkgnd type */
    int		nback_type;						/* nb of params */
    double	back_val[2];					/* user-def. bkg */
    int		nback_val;						/* nb of params */
    int		backsize[2];					/* bkgnd mesh size */
    int		nbacksize;						/* nb of params */
    int		backfsize[2];					/* bkgnd filt. size */
    int		nbackfsize;						/* nb of params */
    double	backfthresh;					/* bkgnd fil. thresh */
    pbacktype			pback_type;	/* phot. bkgnd type */
    int		pback_size;				/* rect. ann. width */
  /*----- memory */
    int		clean_stacksize;			/* size of buffer */
    int		mem_pixstack;				/* pixel stack size */
    int		mem_bufsize;				/* strip height */
  /*----- catalog output */
    char	param_name[MAXCHAR];			/* param. filename */
  /*----- miscellaneous */
    int		pipe_flag;					/* allow piping ? */
    verbosetype      	verbose_type;	/* display type */
    int		xml_flag;				/* Write XML file? */
    char	xml_name[MAXCHAR];			/* XML file name */
    char	xsl_name[MAXCHAR];			/* XSL file name (or URL) */
    char	sdate_start[12];			/* SCAMP start date */
    char	stime_start[12];			/* SCAMP start time */
    char	sdate_end[12];				/* SCAMP end date */
    char	stime_end[12];				/* SCAMP end time */
    double	time_diff;					/* Execution time */
  /*----- CHECK-images */
    int		check_flag;					/* CHECK-image flag */
    checkenum    	check_type[MAXCHECK];		       	/* check-image types */
    int		ncheck_type;				/* nb of params */
    char		*(check_name[MAXCHECK]);       		/* check-image names */
    int		ncheck_name;				/* nb of params */
    struct structcheck	*(check[MAXCHECK]);		/* check-image ptrs */
  /*----- ASSOCiation */
    int		assoc_flag;					/* ASSOCiation flag */
    char		assoc_name[MAXCHAR];	/* ASSOC-list name */
    int		assoc_param[3];				/* ASSOC param cols */
    int		nassoc_param;				/* nb of params */
    int		assoc_data[MAXNASSOC];		/* ASSOC data cols */
    int		nassoc_data;				/* nb of params */

    assoctype	assoc_type;				/* type of assoc. */
    assocselectype  assocselec_type;	/* type of assoc. */
    double	assoc_radius;				/* ASSOC range */
    int		assoc_size;					/* nb of parameters */
    char		retina_name[MAXCHAR];	/* retina filename */
    int		vignetsize[2];				/* vignet size */
    int		vigshiftsize[2];			/* shift-vignet size */
    int		cleanmargin;				/* CLEANing margin */
    char		som_name[MAXCHAR];		/* SOM filename */
    int		somfit_flag;				/* ASSOCiation flag */
    int		somfit_vectorsize;			/* SOMfit vec. size */
  /*----- masking */
    enum 	{MASK_NONE, MASK_BLANK, MASK_CORRECT}	mask_type;	/* type of masking */
    int		blank_flag;					/* BLANKing flag */
    double	weight_thresh[2];      		/* weight threshlds */
    int		nweight_thresh;				/* nb of params */
  /*----- interpolation */
    interptype	interp_type[2];			/* interpolat. type */
    int		ninterp_type;				/* nb of params */
    int		interp_xtimeout[2];   		/* interp. x timeout */
    int		ninterp_xtimeout;       	        /* nb of params */
    int		interp_ytimeout[2];   			/* interp. y timeout */
    int		ninterp_ytimeout;       		/* nb of params */
  /*----- astrometry */
    int		world_flag;				/* WORLD required */
    char		coosys[16];				/* VOTable coord.sys */
    double	epoch;					/* VOTable epoch */
  /*----- growth curve */
    int		growth_flag;				/* gr. curve needed */
    int		flux_growthsize;       			/* number of elem. */
    int		mag_growthsize;       			/* number of elem. */
    int		flux_radiussize;       			/* number of elem. */
    double	growth_step;				/* step size (pix) */
    double	flux_frac[MAXNAPER];			/* for FLUX_RADIUS */
    int		nflux_frac;       			/* number of elem. */
  /*----- PSF-fitting */
    int		psf_flag;				/* PSF-fit needed */
    int		dpsf_flag;				/* dual image PSF-fit */
    char		*(psf_name[2]);				/* PSF filename */
    int		npsf_name;				/* nb of params */
    int		psf_npsfmax;				/* Max # of PSFs */
    psfdisplaytype psfdisplay_type;			/* PSF display type */
    int		psf_xsize,psf_ysize;			/* nb of params */
    int		psf_xwsize,psf_ywsize;			/* nb of params */
    int		psf_alphassize,psf_deltassize;		/* nb of params */
    int		psf_alpha2000size,psf_delta2000size;	/* nb of params */
    int		psf_alpha1950size,psf_delta1950size;	/* nb of params */
    int		psf_fluxsize;				/* nb of params */
    int		psf_fluxerrsize;			/* nb of params */
    int		psf_magsize;				/* nb of params */
    int		psf_magerrsize;				/* nb of params */
    int		pc_flag;				/* PC-fit needed */
    int		pc_vectorsize;				/* nb of params */
    int		prof_flag;				/* Profile-fitting */
    int		pattern_flag;				/* Pattern-fitting */
  /*----- Profile-fitting */
    int		prof_vectorsize;			/* nb of params */
    int		prof_errvectorsize;			/* nb of params */
    int		prof_disk_patternvectorsize;		/* nb of params */
    int		prof_disk_patternncomp;			/* nb of params */
    int		prof_disk_patternmodvectorsize;		/* nb of params */
    int		prof_disk_patternmodncomp;		/* nb of params */
    int		prof_disk_patternargvectorsize;		/* nb of params */
    int		prof_disk_patternargncomp;		/* nb of params */
  /*----- Pattern-fitting */
    pattypenum	pattern_type;				/* Disk pattern type */
  /*----- customize */
    int		fitsunsigned_flag;			/* Force unsign FITS */
    int		next;			     /* Number of extensions in file */
  /* Multithreading */
    int		nthreads;			/* Number of active threads */
    } prefstruct;

/*-------------------------------- catalog  ---------------------------------*/

typedef struct	__attribute__((aligned(16)))
  {
  int		ndetect;				/* nb of detections */
  int		ntotal;					/* Total object nb */
  int		nparam;					/* Nb of parameters */
/*----- Misc. strings defining the extraction */
  char		prefs_name[MAXCHAR];			/* Prefs filename*/
  char		image_name[MAXCHAR];			/* image filename*/
  char		psf_name[MAXCHAR];			/* PSF filename*/
  char		nnw_name[MAXCHAR];			/* NNW name */
  char		filter_name[MAXCHAR];			/* Filter name */
  char		soft_name[MAXCHAR];			/* Sextractor version*/
/*----- time */
  char		ext_date[16],ext_time[16];		/* date and time */
  double	ext_elapsed;				/* processing time */
/*----- MEF */
  int		currext;				/* current extension */
  int		next;					/* Nb of extensions */
  }		sexcatstruct;


  /*------------------------------- structures --------------------------------*/
  typedef struct  __attribute__((aligned(16)))
  	{
  	int	layersnb;
  	int	nn[LAYERS];
  	double	inbias[NEURONS];
  	double	inscale[NEURONS];
  	double	outbias[NEURONS];
  	double	outscale[NEURONS];
  	double	ni[NEURONS];
  	double	no[NEURONS];
  	double	n[LAYERS][NEURONS];
  	double	w[CONNEX][NEURONS][NEURONS];
  	double	b[CONNEX][NEURONS];
  	}	brainstruct;

#endif

/*
 * cudatypes.h
 *
 *  Created on: 10 Jan, 2013
 *      Author: zhao
 */

#ifndef CUDATYPES_H_
#define CUDATYPES_H_

#include "../util5/cudpp.h"
#include "../define.h"

//#define DEBUG_CUDA

#define MAX_THREADS_PER_BLK		256
#define SQUARE_BLK_WIDTH 		16
#define SQUARE_BLK_HEIGHT		16

//#define LARGER_BLK_WIDTH		32
//#define LARGER_BLK_HEIGHT		32

#define	NSONMAX					1024	/* max. number per level */
#define	NBRANCH					16		/* starting number per branch */
#define	MAXPICSIZE				1048576	/* max. image size */

#define	CLEAN_ZONE				10.0	/* zone (in sigma) to */

/*----------------------------object flags---------------------------*/
#define		OBJ_CROWDED		0x0001
#define		OBJ_MERGED		0x0002
#define		OBJ_SATUR		0x0004
#define		OBJ_TRUNC		0x0008
#define		OBJ_APERT_PB	0x0010
#define		OBJ_ISO_PB		0x0020
#define		OBJ_DOVERFLOW	0x0040
#define		OBJ_OVERFLOW	0x0080

/*------------------------background estimation-----------------------*/
#define	BACK_MINGOODFRAC		0.5		/* min frac with good weights*/
#define	QUANTIF_NSIGMA			5		/* histogram limits */
#define	QUANTIF_NMAXLEVELS		4096	/* max nb of quantif. levels */
#define	QUANTIF_AMIN			4		/* min nb of "mode pixels" */
#define	EPS						(1e-4)	/* a small number */

/*--------------------------- Internal constants ----------------------------*/

#define	BIG				1e+30	/* a huge number */
#define	LESSBIG			1e+25	/* a somewhat smaller number */
#define	DATA_BUFSIZE	262144	/* data buffer size */
#define	MARGIN_SCALE	2.0		/* Margin / object height */
#define	MAXCHAR			512		/* max. number of characters */
#define	MAXCHARL		16384	/* max.nb of chars in strlist*/
#define	MAXDEBAREA		3		/* max. area for deblending */
#define	MAXFLAG			4		/* max. # of FLAG-images */
#define	MAXIMAGE		2		/* max. # of input images */
#define	MAXNAPER		32		/* max. number of apertures */
#define	MAXNASSOC		32		/* max. number of assoc. */
#define	MAXPICSIZE		1048576	/* max. image size */
#define	OUTPUT			stderr	/* where all msgs are sent */
#define PSF_NPSFMAX		9		/* Max number of fitted PSFs */
#define EXT_MINAREA		5	 	/* min area in pix to form a object */

#define	DEG			(PI/180.0)	/* 1 deg in radians */
#ifndef PI
#define	PI			3.1415926535898	/* never met before? */
#endif

#define		ANALYSE_FAST		0
#define		ANALYSE_FULL		1
#define		ANALYSE_ROBUST		2

typedef	unsigned int	FLAGTYPE;	/* Flag type */
typedef	unsigned char	BYTE;		/* a byte */
typedef	float		PIXTYPE;		/* Pixel type */

typedef struct
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
  }	objstruct_cuda;

/*---------------data structure to background result-----------*/
extern float 		*d_mean;
extern float 		*d_sigma;


//corresponds to the cdscan field in sextractor
extern float		*d_cdPixArray; //cdvalue field

//corresponds to the the scan field in sextractor
extern float		*d_pixelArray; //dvalue field
extern int			*d_labelArray;
extern int			*d_equivArray;
extern unsigned int	*d_indexArray;


extern unsigned int	*d_compactMask;
extern size_t		*d_numPixAboveThresh;

////////////////temp variables, allocated in the first detection pass
extern unsigned int	*d_segmentMask;
extern unsigned int	*d_compactedIndexArray;
extern unsigned int	*d_compactedLabelArray;
extern unsigned int	*d_compactedParentLabel;
extern float		*d_compactedcdPixelArray;
extern float		*d_compactedPixelArray;
extern unsigned int	*d_pixelCountMask;
extern unsigned int	*d_pixelCountSegment;

extern float		*d_fdpeakSegment;
extern float		*d_fdfluxSegment;
extern unsigned int *d_prunedSegmentMask;

//The output of the first scan.
extern unsigned int	*d_pixelCountArray;
extern float		*d_fdpeakArray;
extern float		*d_fdfluxArray;
extern float		*d_dthreshArray;
extern unsigned int	*d_ok;

extern unsigned int	*d_finalPixelIndexArray;
extern unsigned int	*d_finalLabelArray;
extern unsigned int	*d_finalObjIndexArray;

/*-----------data structure to store objects after cutting-----------*/
extern unsigned int	*d_cuttedObjLevelArray;
extern unsigned int	*d_cuttedObjLabelArray;
extern unsigned int	*d_cuttedObjIndexArray;
extern unsigned int	*d_cuttedObjFlagArray;
extern unsigned int	*d_cuttedPixCountArray;
extern unsigned int	*d_cuttedRootlabelArray;
extern float		*d_cuttedDthreshArray;

/*-----------data structure to store objects attributes-----------*/
extern unsigned int	*d_index;
extern unsigned int	*d_flag;
extern unsigned int	*d_fdnpix;
extern float		*d_dthresh;

extern unsigned int *d_xmin;
extern unsigned int *d_xmax;
extern unsigned int *d_ymin;
extern unsigned int *d_ymax;
extern unsigned int *d_dnpix;
extern unsigned int *d_npix;
extern unsigned int	*d_peakx;
extern unsigned int	*d_peaky;

extern float 		*d_bkg;
extern float 		*d_dbkg;
extern float		*d_sigbkg;
extern float 		*d_a;
extern float 		*d_b;
extern float 		*d_cxx;
extern float 		*d_cxy;
extern float 		*d_cyy;
extern float 		*d_theta;
extern float 		*d_abcor;
extern float		*d_peak;
extern float 		*d_dpeak;
extern float 		*d_fdpeak;
extern float		*d_flux;
extern float 		*d_dflux;
extern float 		*d_fdflux;
extern float		*d_fluxerr;
extern float		*d_thresh;
extern float		*d_mthresh;
extern float		*d_fwhm;

extern double 		*d_mx;
extern double 		*d_my;
extern double 		*d_mx2;
extern double 		*d_my2;
extern double 		*d_mxy;
extern double		*d_poserr_mx2;
extern double		*d_poserr_my2;
extern double		*d_poserr_mxy;
extern char  		*d_singuflag;

extern int			*d_iso[NISO];

//////////////////////////
//used in endobject phase
extern double	*d_mxw;
extern double	*d_myw;
extern double	*d_alphas;
extern double	*d_deltas;
extern double	*d_alpha2000;
extern double	*d_delta2000;

extern double	*d_posx;
extern double	*d_posy;
extern float	*d_elong;
extern float	*d_ellip;
extern float	*d_polar;
extern float	*d_sprob;
extern float	*d_sposx;
extern float	*d_sposy;
//extern float	*d_poserr_a;
//extern float	*d_poserr_b;
//extern float	*d_poserr_theta;
//extern float	*d_poserr_cxx;
//extern float	*d_poserr_cyy;
//extern float	*d_poserr_cxy;
extern float	*d_flux_iso;
extern float	*d_fluxerr_iso;
extern float	*d_flux_isocor;
extern float	*d_fluxerr_isocor;
extern float	*d_mag_iso;
extern float	*d_magerr_iso;
extern float	*d_kronfactor;
extern float	*d_flux_auto;
extern float	*d_fluxerr_auto;
extern float	*d_mag_auto;
extern float	*d_magerr_auto;

extern float	*d_flux_aper[NAPER];
extern float	*d_fluxerr_aper[NAPER];
extern float	*d_mag_aper[NAPER];
extern float	*d_magerr_aper[NAPER];

/////////////////////////

extern CUDPPHandle 	theCudpp;

extern CUDPPConfiguration config;
extern CUDPPHandle 	scanplan;
extern CUDPPHandle	reduceplan;
extern CUDPPHandle	segscanplan;
extern CUDPPHandle	compactplan;
extern CUDPPHandle	sortplan;


extern int 		width;
extern int 		height;
extern float 	thresh;
extern float 	global_dthresh;

extern cudaEvent_t start_t, stop_t;
extern cudaEvent_t start, stop;

extern unsigned int *d_masterIndex;

#endif /* CUDATYPES_H_ */

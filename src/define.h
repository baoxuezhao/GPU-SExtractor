 /*
 				define.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global definitions.
*
*	Last modify:	31/03/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
#ifndef __DEFINE_H__
#define __DEFINE_H__

/* Check if we are using a configure script here */
#ifndef HAVE_CONFIG_H
#define		VERSION		"2.x"
#define		DATE		"2009-03-31"
#define		THREADS_NMAX	16		/* max. number of threads */
#endif

/*------------------------ what, who, when and where ------------------------*/

#define		BANNER		"SExtractor"
#define		MYVERSION	VERSION
#define		EXECUTABLE	"sex"
#define		COPYRIGHT	"Emmanuel BERTIN <bertin@iap.fr>"
#define		WEBSITE		"http://astromatic.iap.fr/software/sextractor"
#define		INSTITUTE	"IAP  http://www.iap.fr"

/*--------------------------- Internal constants ----------------------------*/

#define	BIG				1e+30		/* a huge number */
#define	LESSBIG			1e+25		/* a somewhat smaller number */
#define	DATA_BUFSIZE	262144		/* data buffer size */
#define	MARGIN_SCALE	2.0			/* Margin / object height */
#define	MAXCHAR			512			/* max. number of characters */
#define	MAXCHARL		16384		/* max.nb of chars in strlist*/
#define	MAXDEBAREA		3			/* max. area for deblending */
#define	MAXFLAG			4			/* max. # of FLAG-images */
#define	MAXIMAGE		2			/* max. # of input images */
#define	MAXNAPER		32			/* max. number of apertures */
#define	MAXNASSOC		32			/* max. number of assoc. */
#define	MAXPICSIZE		1048576		/* max. image size */
#define	NISO			8			/* number of isophotes */
#define	NAPER			1		/* number of apertures */
#define	OUTPUT			stderr		/* where all msgs are sent */
#define PSF_NPSFMAX		9			/* Max number of fitted PSFs */

#define	DEG				(PI/180.0)	/* 1 deg in radians */
#ifndef PI
#define	PI				3.1415926535898	/* never met before? */
#endif

/* NOTES:
 *
 *One must have:	BIG < the biggest element a float can store
 *			DATA_BUFSIZE >= 2880 with DATA_BUFSIZE%8 = 0
 *			MAXCHAR >= 16
 *			1 <= MAXCHECK <= MAXLIST (see prefs.h)
 *			1 <= MAXDEBAREA (see prefs.c & extract.c)
 *			1 <= MAXFLAG <= MAXLIST (see prefs.h)
 *			1 <= MAXIMAGE <= MAXLIST (see prefs.h)
 *			1 <= MAXNAPER <= MAXLIST (see prefs.h)
 *			1 <= MAXNASSOC <= MAXLIST (see prefs.h)
 *			MAXPICSIZE > size of any image!!
 *			NISO = 8 (otherwise need to change prefs.h)
 *			1 <= PSF_NPSFMAX
*/

/*--------------------------- Neural Network parameters ---------------------*/
#define		LAYERS		3	/* max. number of hidden+i/o layers */
#define		CONNEX		LAYERS-1
#define		NEURONS		10	/* maximum number of neurons/layer */


/*---- Set defines according to machine's specificities and customizing -----*/

#if _LARGEFILE_SOURCE
#define	FSEEKO	fseeko
#define	FTELLO	ftello
#else
#define	FSEEKO	fseek
#define	FTELLO	ftell
#endif
/*--------------------- in case of missing constants ------------------------*/

#ifndef		SEEK_SET
#define		SEEK_SET	0
#endif
#ifndef		SEEK_CUR
#define		SEEK_CUR	1
#endif

#ifndef	EXIT_SUCCESS
#define	EXIT_SUCCESS		0
#endif
#ifndef	EXIT_FAILURE
#define	EXIT_FAILURE		-1
#endif

 #define	MAXLIST		32		/* max. nb of list members */
/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

#if defined (__CUDACC__)
#define ALIGN16       __align__(16)
#elif defined (__GCCXML__)
#define ALIGN16
#elif defined (_MSC_VER)
#define ALIGN16         __declspec(align(16))
#elif defined (__GNUC__)
#define ALIGN16         __attribute__((__aligned__(16)))
#else
#error compiler unsupport for alignment
#endif

/*-------------------------------- flags ------------------------------------*/

#define		OBJ_CROWDED	0x0001
#define		OBJ_MERGED	0x0002
#define		OBJ_SATUR	0x0004
#define		OBJ_TRUNC	0x0008
#define		OBJ_APERT_PB	0x0010
#define		OBJ_ISO_PB	0x0020
#define		OBJ_DOVERFLOW	0x0040
#define		OBJ_OVERFLOW	0x0080

/*----------------------------- weight flags --------------------------------*/

#define		OBJ_WEIGHTZERO	0x0001
#define		OBJ_DWEIGHTZERO 0x0002

/*---------------------------- preanalyse flags -----------------------------*/

#define		ANALYSE_FAST		0
#define		ANALYSE_FULL		1
#define		ANALYSE_ROBUST		2

/*----------------------------- FITS WCS --------------------------*/

#define		NAXIS	2		/* Max number of FITS axes */

#define		DEG	(PI/180.0)	/* 1 deg in radians */
#define		ARCMIN	(DEG/60.0)	/* 1 arcsec in radians */
#define		ARCSEC	(DEG/3600.0)	/* 1 arcsec in radians */
#define		MAS	(ARCSEC/1000.0)	/* 1 mas in radians */
#define		MJD2000	51544.50000	/* Modified Julian date for J2000.0 */
#define		MJD1950	33281.92346	/* Modified Julian date for B1950.0 */
#define		JU2TROP	1.0000214	/* 1 Julian century in tropical units*/
#define		WCS_NOCOORD	1e31	/* Code for non-existing coordinates */

#define		WCS_NGRIDPOINTS	12	/* Number of WCS grid points / axis */
#define		WCS_NGRIDPOINTS2	(WCS_NGRIDPOINTS*WCS_NGRIDPOINTS)
#define		WCS_INVMAXDEG	9	/* Maximum inversion polynom degree */
#define		WCS_INVACCURACY	0.04	/* Maximum inversion error (pixels) */
#define		WCS_NRANGEPOINTS 32	/* Number of WCS range points / axis */

#define 	WCSSET 137
#define 	LINSET 137

#define 	WCSTRIG_TOL 1e-10
/*---------------------------- photometry -----------------------------*/

#define	APER_OVERSAMP	5	/* oversampling in each dimension (MAG_APER) */
#define	KRON_NSIG	3*MARGIN_SCALE	/* MAG_AUTO analysis range (number */
					/* of sigma) */
#define	PETRO_NSIG	3*MARGIN_SCALE	/* MAG_PETRO analysis range (number */
					/* of sigma) */
#define	CROWD_THRESHOLD	0.1	/* The OBJ_CROWDED flag is set if photometric*/
				/* contamination may exceed this fraction of */
				/* flux */

/* NOTES:
One must have:	APER_OVERSAMP >= 1
		KRON_NSIG > 0.0
		PETRO_NSIG > 0.0
		CROWD_THRESHOLD >= 0
*/

/*--------------------------- Neural Network parameters ---------------------*/
#define		LAYERS		3	/* max. number of hidden+i/o layers */
#define		CONNEX		LAYERS-1
#define		NEURONS		10	/* maximum number of neurons/layer */


/*------------------- a few definitions to read FITS parameters ------------*/

#define	FBSIZE	2880L	/* size (in bytes) of one FITS block */

#define	FITSTOF(k, def) \
			(st[0]=0,((point = fitsnfind(buf, k, n))? \
				 fixexponent(point), \
				atof(strncat(st, &point[10], 70)) \
				:(def)))

#define	FITSTOI(k, def) \
			(st[0]=0,(point = fitsnfind(buf, k, n))? \
				 atoi(strncat(st, &point[10], 70)) \
				:(def))

#define	FITSTOS(k, str, def) \
                { if (fitsread(buf,k,str,H_STRING,T_STRING)!= RETURN_OK) \
                    strcpy(str, (def)); \
                }

/*------------------------------- Other Macros -----------------------------*/

#define	DEXP(x)	exp(2.30258509299*(x))	/* 10^x */

#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define	QFSEEK(afile, offset, pos, fname) \
		if (FSEEKO(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#define	QFTELL(afile, pos, fname) \
		if ((pos=FTELLO(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QFREE(ptr) \
		{free(ptr); \
		ptr = NULL;}

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ)))) \
		   error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
                    error(EXIT_FAILURE, "Not enough memory for ", \
                        #ptrout " (" #nel " elements) !"); \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ));};}

#define	RINT(x)	(int)(floor(x+0.5))

#define	PIX(pic, x, y)	pic->strip[(((int)y)%pic->stripheight) \
				*pic->width +(int)x]

#define	NPRINTF		if (prefs.verbose_type == NORM \
				|| prefs.verbose_type==WARN) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type==NORM \
				|| prefs.verbose_type==WARN) \
				fprintf(w, "\33[1M> %s\n\33[1A",x); \
			else if (prefs.verbose_type == FULL) \
				fprintf(w, "%s.\n", x);}

#define	QPRINTF		if (prefs.verbose_type != QUIET)	fprintf

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf

#define	QWARNING       	if (prefs.verbose_type==WARN \
				|| prefs.verbose_type==FULL)	warning

#define FLAG(x)		(*((char *)&flag##x))

#define VECFLAG(x)	(*((char *)flag##x))

#endif

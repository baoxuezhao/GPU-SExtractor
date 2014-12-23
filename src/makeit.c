/*
 				makeit.c

 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	Part of:	SExtractor
 *
 *	Author:		E.BERTIN, IAP & Leiden observatory
 *
 *	Contents:	main program.
 *
 *	Last modify:	20/11/2008
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"assoc.h"
#include	"back.h"
#include	"check.h"
#include	"fft.h"
#include	"field.h"
#include	"filter.h"
#include	"growth.h"
#include	"interpolate.h"
#include	"pattern.h"
#include	"psf.h"
#include	"profit.h"
#include	"som.h"
#include	"weight.h"
#include	"xml.h"

#include	"cuda/cudainit.h"
#include 	"cuda/cudaback.h"

static int		selectext(char *filename);
time_t			thetimet, thetimet2;
extern profitstruct	*theprofit;
extern char		profname[][32];

/******************************** makeit *************************************/
/*
Manage the whole stuff.
 */
void	makeit()
{

	checkstruct		*check;
	picstruct		*field, *pffield[MAXFLAG];
	catstruct		*imacat;
	tabstruct		*imatab;
	patternstruct	*pattern;
	static time_t        thetime1, thetime2;
	struct tm		*tm;
	int			nflag[MAXFLAG],
	i, nok, ntab, next, ntabmax, forcextflag,
	nima0,nima1, nweight0,nweight1, npat;

	/* Install error logging */
	//error_installfunc(write_error);

	/* Processing start date and time */
	thetimet = time(NULL);
	tm = localtime(&thetimet);
	sprintf(prefs.sdate_start,"%04d-%02d-%02d",
			tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
	sprintf(prefs.stime_start,"%02d:%02d:%02d",
			tm->tm_hour, tm->tm_min, tm->tm_sec);

	clock_t _start = clock();

	NFPRINTF(OUTPUT, "");
	QPRINTF(OUTPUT, "----- %s %s started on %s at %s with %d thread%s\n\n",
			BANNER,
			MYVERSION,
			prefs.sdate_start,
			prefs.stime_start,
			prefs.nthreads,
			prefs.nthreads>1? "s":"");

	/* Initialize globals variables */
	initglob();

	NFPRINTF(OUTPUT, "Setting catalog parameters");
	readcatparams(prefs.param_name);
	useprefs();			/* update things accor. to prefs parameters */

	if (prefs.filter_flag)
	{
		NFPRINTF(OUTPUT, "Reading detection filter");
		getfilter(prefs.filter_name);	/* get the detection filter */
	}

	if (FLAG(obj2.sprob))
	{
		NFPRINTF(OUTPUT, "Initializing Neural Network");
		neurinit();
		NFPRINTF(OUTPUT, "Reading Neural Network Weights");
		getnnw();
	}

	/* Allocate memory for multidimensional catalog parameter arrays */
	alloccatparams();

	useprefs();

	/* Check if a specific extension should be loaded */
	if ((nima0=selectext(prefs.image_name[0])) != RETURN_ERROR)
	{
		forcextflag = 1;
		ntabmax = next = 1;
	}
	else
		forcextflag = 0;

	if (!(imacat = read_cat(prefs.image_name[0])))
		error(EXIT_FAILURE, "*Error*: cannot open ", prefs.image_name[0]);
	close_cat(imacat);
	imatab = imacat->tab;

	if (!forcextflag)
	{
		ntabmax = imacat->ntab;
		/*-- Compute the number of valid input extensions */
		next = 0;
		for (ntab = 0 ; ntab<imacat->ntab; ntab++, imatab = imatab->nexttab)
		{
			/*---- Check for the next valid image extension */
			if ((imatab->naxis < 2)
					|| !strncmp(imatab->xtension, "BINTABLE", 8)
					|| !strncmp(imatab->xtension, "ASCTABLE", 8))
				continue;
			next++;
		}
	}

	/* Do the same for other data (but do not force single extension mode) */
	nima1 = selectext(prefs.image_name[1]);
	nweight0 = selectext(prefs.wimage_name[0]);
	nweight1 = selectext(prefs.wimage_name[1]);
	for (i=0; i<prefs.nfimage_name; i++)
		nflag[i] = selectext(prefs.fimage_name[i]);

	thecat.next = next;

	NFPRINTF(OUTPUT, "Initializing catalog");
	initcat();

	/* Go through all images */
	nok = -1;
	for (ntab = 0 ; ntab<ntabmax; ntab++, imatab = imatab->nexttab)
	{
		/*--  Check for the next valid image extension */
		if (!forcextflag && ((imatab->naxis < 2)
				|| !strncmp(imatab->xtension, "BINTABLE", 8)
				|| !strncmp(imatab->xtension, "ASCTABLE", 8)))
			continue;
		nok++;

		/*-- Initial time measurement*/
		time(&thetime1);
		thecat.currext = nok+1;

		field = newfield(prefs.image_name[0], DETECT_FIELD | MEASURE_FIELD,
				nima0<0? nok:nima0);

		/* read all the data into the image buffer*/
		field->strip = (float*)malloc(field->width*field->height*sizeof(PIXTYPE));
		field->ymax = field->height;

		read_body(field->tab, field->strip, field->width*field->height);

		QPRINTF(OUTPUT, "Initialize the GPU device: ");


		init_device(field->width, field->height, field->strip);


		/*-- Compute background maps for `standard' fields */
		QPRINTF(OUTPUT, "Detection+Measurement image: ");

		makeback(field); //0.28s

		QPRINTF(OUTPUT,  "(M+D) Background: %-10g RMS: %-10g / Threshold: %-10g \n",
				field->backmean, field->backsig, (field->flags & DETECT_FIELD)?
						field->dthresh: field->thresh);

		thefield1 = thefield2 = *field;

		reinitcat(field);

		/*-- Start the extraction pipeline */
		NFPRINTF(OUTPUT, "Scanning image");
		scanimage(field, pffield, prefs.nimaflag);

		/*-- Final time measurements*/
		if (time(&thetime2)!=-1)
		{
			if (!strftime(thecat.ext_date, 12, "%d/%m/%Y", localtime(&thetime2)))
				error(EXIT_FAILURE, "*Internal Error*: Date string too long ","");
			if (!strftime(thecat.ext_time, 10, "%H:%M:%S", localtime(&thetime2)))
				error(EXIT_FAILURE, "*Internal Error*: Time/date string too long ","");
			thecat.ext_elapsed = difftime(thetime2, thetime1);
		}

		reendcat();

		/*-- Close ASSOC routines */
		//end_assoc(field);

		endfield(field);

		QPRINTF(OUTPUT, "Objects: detected %-8d / sextracted %-8d               \n",
				thecat.ndetect, thecat.ntotal);
	}

	if (nok<0)
		error(EXIT_FAILURE, "Not enough valid FITS image extensions in ",
				prefs.image_name[0]);
	free_cat(&imacat, 1);

	NFPRINTF(OUTPUT, "Closing files");

	if (prefs.filter_flag)
		endfilter();

	if (FLAG(obj2.sprob))
		neurclose();

	/* Processing end date and time */
	thetimet2 = time(NULL);
	tm = localtime(&thetimet2);
	sprintf(prefs.sdate_end,"%04d-%02d-%02d",
			tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
	sprintf(prefs.stime_end,"%02d:%02d:%02d",
			tm->tm_hour, tm->tm_min, tm->tm_sec);
	prefs.time_diff = difftime(thetimet2, thetimet);

	endcat((char *)NULL);

	clock_t _end = clock();
	printf("time consumed by the SExtractor (including output) is %f (sec)\n", (_end-_start+0.0)/CLOCKS_PER_SEC);

	return;
}


/******************************** initglob ***********************************/
/*
Initialize a few global variables
 */
void	initglob()
{
	int	i;

	for (i=0; i<37; i++)
	{
		ctg[i] = cos(i*PI/18);
		stg[i] = sin(i*PI/18);
	}


	return;
}


/****** selectext ************************************************************
PROTO	int selectext(char *filename)
PURPOSE	Return the user-selected extension number [%d] from the file name.
INPUT	Filename character string.
OUTPUT	Extension number, or RETURN_ERROR if nos extension specified.
NOTES	The bracket and its extension number are removed from the filename if
	found.
AUTHOR  E. Bertin (IAP)
VERSION 08/10/2007
 ***/
static int	selectext(char *filename)
{
	char	*bracl,*bracr;
	int	next;

	if (filename && (bracl=strrchr(filename, '[')))
	{
		*bracl = '\0';
		if ((bracr=strrchr(bracl+1, ']')))
			*bracr = '\0';
		next = strtol(bracl+1, NULL, 0);
		return next;
	}

	return RETURN_ERROR;
}


/****** write_error ********************************************************
PROTO	int	write_error(char *msg1, char *msg2)
PURPOSE	Manage files in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/07/2006
 ***/
void	write_error(char *msg1, char *msg2)
{
	char			error[MAXCHAR];

	sprintf(error, "%s%s", msg1,msg2);
	if (prefs.xml_flag)
		write_xmlerror(prefs.xml_name, error);

	/* Also close existing catalog */
	endcat(error);

	end_xml();

	return;
}


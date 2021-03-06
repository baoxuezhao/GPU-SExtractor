/*
 readimage.c

 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	Part of:	SExtractor
 *
 *	Author:		E.BERTIN (IAP)
 *
 *	Contents:	functions for input of image data.
 *
 *	Last modify:	11/10/2007
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

#include	"wcs/wcs.h"
#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"check.h"
#include	"field.h"
#include	"fits/fitscat.h"
#include	"fitswcs.h"
#include	"interpolate.h"
#include	"back.h"
#include	"astrom.h"
#include	"weight.h"
#include    "wcs/tnx.h"

#include	"cuda/cudaback.h"
/******************************* loadstrip ***********************************/
/*
 Load a new strip of pixel data into the buffer.
 */
void *loadstrip(picstruct *field)

{
	tabstruct *tab;
	checkstruct *check;
	int y, w, flags, interpflag;
	PIXTYPE *data, *wdata, *rmsdata;

	tab = field->tab;
	w = field->width;
	flags = field->flags;
	interpflag = 0; //(wfield && wfield->interp_flag);
	wdata = NULL; /* To avoid gcc -Wall warnings */

	if (!field->y)
	{
		/*- First strip */
		int nbpix;

		nbpix = w * field->stripheight;

		if (flags ^ FLAG_FIELD) {

			/*---- Allocate space for the frame-buffer */
			/*
			if (!(field->strip = (PIXTYPE *) malloc(
					field->stripheight * field->width * sizeof(PIXTYPE))))
				error(EXIT_FAILURE,
						"Not enough memory for the image buffer of ",
						field->rfilename); */

			data = field->strip;
			/*---- We assume weight data have been read just before */
			/*
			if (interpflag)
				wdata = wfield->strip;
			if (flags & BACKRMS_FIELD)
				for (y = 0, rmsdata = data; y < field->stripheight;
						y++, rmsdata += w)
					backrmsline(field, y, rmsdata);
			else if (flags & INTERP_FIELD)
				copydata(field, 0, nbpix);
			else *///true

				//read_body(tab, data, nbpix);

			//flag is in DETECT_FIELD | MEASURE_FIELD
			//if (flags & (WEIGHT_FIELD | RMS_FIELD | BACKRMS_FIELD | VAR_FIELD)) //false
			//	weight_to_var(field, data, nbpix);

			if ((flags & MEASURE_FIELD)
					&& (check = prefs.check[CHECK_IDENTICAL]))	//false
				writecheck(check, data, nbpix);

			for (y = 0; y < field->stripheight; y++, data += w) {
				/*------ This is the only place where one can pick-up safely the current bkg */
				if (flags & (MEASURE_FIELD | DETECT_FIELD))
					subbackline(field, y, data);

				/*------ Go to interpolation process */
				/*
				if (interpflag) {
					interpolate(field, wfield, data, wdata);
					wdata += w;
				}*/
				/*------ Check-image stuff */
				if (prefs.check_flag) {		//false

					if (flags & MEASURE_FIELD) {
						if ((check = prefs.check[CHECK_BACKGROUND]))
							writecheck(check, field->backline, w);
						if ((check = prefs.check[CHECK_SUBTRACTED]))
							writecheck(check, data, w);
						if ((check = prefs.check[CHECK_APERTURES]))
							writecheck(check, data, w);
						if ((check = prefs.check[CHECK_SUBPSFPROTOS]))
							writecheck(check, data, w);
						if ((check = prefs.check[CHECK_SUBPCPROTOS]))
							writecheck(check, data, w);
						if ((check = prefs.check[CHECK_SUBPROFILES]))
							writecheck(check, data, w);
					}
					if ((flags & DETECT_FIELD) && (check =
							prefs.check[CHECK_BACKRMS])) {  //false

						backrmsline(field, y, (PIXTYPE *) check->pix);
						writecheck(check, (PIXTYPE*)(check->pix), w);
					}
				}
			}
		}
		else {
			if (!(field->fstrip = (FLAGTYPE *) malloc(
					field->stripheight * field->width * sizeof(FLAGTYPE))))
				error(EXIT_FAILURE, "Not enough memory for the flag buffer of ",
						field->rfilename);
			read_ibody(field->tab, field->fstrip, nbpix);
		}

		field->ymax = field->stripheight;

		if (field->ymax < field->height)
			field->stripysclim = field->stripheight - field->stripmargin;
	} else {

		/*- other strips */
		if (flags ^ FLAG_FIELD) {
			data = field->strip + field->stripylim * w;
			/*---- We assume weight data have been read just before */
			//if (interpflag)
			//	wdata = wfield->strip + field->stripylim * w;

			/*---- copy to Check-image the "oldest" line before it is replaced */
			if ((flags & MEASURE_FIELD) && (check =
					prefs.check[CHECK_SUBOBJECTS]))
				writecheck(check, data, w);

			if (flags & BACKRMS_FIELD)
				backrmsline(field, field->ymax, data);
			else if (flags & INTERP_FIELD)
				copydata(field, field->stripylim * w, w);
			else
				read_body(tab, data, w);

			//if (flags & (WEIGHT_FIELD | RMS_FIELD | BACKRMS_FIELD | VAR_FIELD))
			//	weight_to_var(field, data, w);

			if ((flags & MEASURE_FIELD)
					&& (check = prefs.check[CHECK_IDENTICAL]))
				writecheck(check, data, w);
			/*---- Interpolate and subtract the background at current line */
			if (flags & (MEASURE_FIELD | DETECT_FIELD))
				subbackline(field, field->ymax, data);

			//if (interpflag)
			//	interpolate(field, wfield, data, wdata);
			/*---- Check-image stuff */
			if (prefs.check_flag) {
				if (flags & MEASURE_FIELD) {
					if ((check = prefs.check[CHECK_BACKGROUND]))
						writecheck(check, field->backline, w);
					if ((check = prefs.check[CHECK_SUBTRACTED]))
						writecheck(check, data, w);
					if ((check = prefs.check[CHECK_APERTURES]))
						writecheck(check, data, w);
					if ((check = prefs.check[CHECK_SUBPSFPROTOS]))
						writecheck(check, data, w);
					if ((check = prefs.check[CHECK_SUBPCPROTOS]))
						writecheck(check, data, w);
					if ((check = prefs.check[CHECK_SUBPROFILES]))
						writecheck(check, data, w);
				}
				if ((flags & DETECT_FIELD) && (check =
						prefs.check[CHECK_BACKRMS])) {
					backrmsline(field, field->ymax, (PIXTYPE *) check->pix);
					writecheck(check, (PIXTYPE*)(check->pix), w);
				}
			}
		} else
			read_ibody(tab, field->fstrip + field->stripylim * w, w);

		field->stripylim = (++field->ymin) % field->stripheight;
		if ((++field->ymax) < field->height)
			field->stripysclim = (++field->stripysclim) % field->stripheight;
	}

	return (flags ^ FLAG_FIELD) ?
			(void *) (field->strip + field->stripy * w) :
			(void *) (field->fstrip + field->stripy * w);
}

/******************************** copydata **********************************/
/*
 Copy image data from one field to the other.
 */
void copydata(picstruct *field, int offset, int size) {
	memcpy(field->strip + offset, field->reffield->strip + offset,
			size * sizeof(PIXTYPE));
	return;
}

/******************************* readimagehead *******************************/
/*
 extract some data from the FITS-file header
 */
void readimagehead(picstruct *field) {
#define FITSREADS(buf, k, str, def) \
                {if (fitsread(buf,k,str, H_STRING,T_STRING) != RETURN_OK) \
                   strcpy(str, (def)); \
                }

	tabstruct *tab;

	tab = field->tab;

	if (tab->naxis < 2)
		error(EXIT_FAILURE, field->filename, " does NOT contain 2D-data!");

	/*---------------------------- Basic keywords ------------------------------*/
	if (tab->bitpix != BP_BYTE && tab->bitpix != BP_SHORT
			&& tab->bitpix != BP_LONG && tab->bitpix != BP_FLOAT
			&& tab->bitpix != BP_DOUBLE)
		error(EXIT_FAILURE, "Sorry, I don't know that kind of data.", "");

	field->width = tab->naxisn[0];
	field->height = tab->naxisn[1];
	field->npix = (KINGSIZE_T) field->width * field->height;
	if (tab->bitsgn && prefs.fitsunsigned_flag)
		tab->bitsgn = 0;

	FITSREADS(tab->headbuf, "OBJECT  ", field->ident, "Unnamed");

	/*----------------------------- Astrometry ---------------------------------*/
	/* Presently, astrometry is done only on the measurement and detect images */
	if (field->flags & (MEASURE_FIELD | DETECT_FIELD))
		field->wcs = read_wcs(tab);

	QFSEEK(field->file, tab->bodypos, SEEK_SET, field->filename);

	return;

#undef FITSREADS

}


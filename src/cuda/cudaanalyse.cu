/*
 * cudaanalyse.cu
 *
 *  Created on: 17 Mar, 2013
 *      Author: zhao
 */
#include <stdio.h>
#include "cudatypes.h"

#include "../types.h"
#include "../globals.h"

//prenalysis before gather up, for level 0 to n.
// level 0 takes almost half the whole time
//thresh is a constant here


__device__ float hmedian1(float *ra, int n)
{
	int		l, j, ir, i;
	float	rra;

	if (n<2)
		return *ra;
	ra--;
	for (l = ((ir=n)>>1)+1;;)
	{
		if (l>1)
			rra = ra[--l];
		else
		{
			rra = ra[ir];
			ra[ir] = ra[1];
			if (--ir == 1)
			{
				ra[1] = rra;
				return n&1? ra[n/2+1] : (float)((ra[n/2]+ra[n/2+1])/2.0);
			}
		}
		for (j = (i=l)<<1; j <= ir;)
		{
			if (j < ir && ra[j] < ra[j+1])
				++j;
			if (rra < ra[j])
			{
				ra[i] = ra[j];
				j += (i=j);
			}
			else
				j = ir + 1;
		}
		ra[i] = rra;
	}

	/* (the 'return' is inside the loop!!) */
}

/****************************** computeisocorflux ****************************/
/*
Compute the (corrected) isophotal flux.
 */
__device__ void  computeisocorflux(
		unsigned int 	*d_fdnpix,
		float			*d_dthresh,
		float 			*d_flux,
		float			*d_fluxerr,
		float			*d_flux_isocor,
		float			*d_fluxerr_isocor,
		float			flagobj2_fluxerr_isocor,
		int tid)

{
	double	ati;

	ati = (d_flux[tid]>0.0)? (d_fdnpix[tid]*d_dthresh[tid]/d_flux[tid]) : 0.0;
	if (ati>1.0)
		ati = 1.0;
	else if (ati<0.0)
		ati = 0.0;
	d_flux_isocor[tid] = d_flux[tid]/(1.0-0.196099*ati-0.751208*ati*ati);

	if (*((char *)&flagobj2_fluxerr_isocor))
	{
		if (d_flux[tid] > 0.0)
		{
			double	dati, sigtv;

			sigtv = d_fluxerr[tid]/(d_flux[tid]*d_flux[tid]);
			dati = d_fdnpix[tid]?ati*sqrt(sigtv+1.0/d_fdnpix[tid]): 0.0;
			dati = 0.196099*dati + 0.751208*2*ati*dati;
			d_fluxerr_isocor[tid] = sqrt(sigtv+dati*dati)*d_flux[tid];
		}
		else
			d_fluxerr_isocor[tid] = sqrt(d_fluxerr[tid]);
	}

	return;
}

/***************************** computeaperflux********************************/
/*
Compute the total flux within a circular aperture.
 */
__device__ void  computeaperflux(

		float		*d_pixelArray,
		int			*d_labelArray, //allocation map
		unsigned int *d_flag,
		float		*d_dbkg,
		float		**d_flux_aper,
		float		**d_fluxerr_aper,
		float		**d_mag_aper,
		float		**d_magerr_aper,
		double		*d_mx,
		double		*d_my,
		prefstruct	*d_prefs,
		int 		width,
		int			height,
		int			fymin,
		int			fymax,
		double		ngamma,
		double		gain,
		float		backsig,
		int			tid,
		int 		i)

{
	float		r2, raper,raper2, rintlim,rintlim2,rextlim2, dx,dx1,dy,dy2,
	offsetx,offsety,scalex,scaley,scale2, locarea;
	double		tv, sigtv, area, pix, var, backnoise2;
	int			x,y, x2,y2, xmin,xmax,ymin,ymax, sx,sy, pflag,corrflag;
	long			pos;
	PIXTYPE		*strip,*stript;

	//if (wfield)
	//	wthresh = wfield->weight_thresh;
	//wstrip = wstript = NULL;
	//mx = obj->mx;
	//my = obj->my;
	//w = field->width;
	//h = field->stripheight;
	//fymin = field->ymin;
	//fymax = field->ymax;
	//ngamma = field->ngamma;
	pflag = (d_prefs->detect_type==d_prefs->PHOTO)? 1:0;
	corrflag = (d_prefs->mask_type==d_prefs->MASK_CORRECT);
	//gainflag = wfield && prefs.weightgain_flag;
	var = backnoise2 = backsig*backsig;
	//gain = field->gain;
	/* Integration radius */
	raper = d_prefs->apert[i]/2.0;
	raper2 = raper*raper;
	/* Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
	rintlim = raper - 0.75;
	rintlim2 = (rintlim>0.0)? rintlim*rintlim: 0.0;
	/* External radius of the oversampled annulus (>r+sqrt(2)/2) */
	rextlim2 = (raper + 0.75)*(raper + 0.75);
	tv = sigtv = area = 0.0;
	scaley = scalex = 1.0/APER_OVERSAMP;
	scale2 = scalex*scaley;
	offsetx = 0.5*(scalex-1.0);
	offsety = 0.5*(scaley-1.0);

	xmin = (int)(d_mx[tid]-raper+0.499999);
	xmax = (int)(d_mx[tid]+raper+1.499999);
	ymin = (int)(d_my[tid]-raper+0.499999);
	ymax = (int)(d_my[tid]+raper+1.499999);

	if (xmin < 0)
	{
		xmin = 0;
		d_flag[tid] |= OBJ_APERT_PB;
	}
	if (xmax > width)
	{
		xmax = width;
		d_flag[tid] |= OBJ_APERT_PB;
	}
	if (ymin < fymin)
	{
		ymin = fymin;
		d_flag[tid] |= OBJ_APERT_PB;
	}
	if (ymax > fymax)
	{
		ymax = fymax;
		d_flag[tid] |= OBJ_APERT_PB;
	}

	strip = d_pixelArray;
	//if (wfield)
	//	wstrip = wfield->strip;
	for (y=ymin; y<ymax; y++)
	{
		stript = strip + (pos = (y%height)*width + xmin);
		//if (wfield)
		//	wstript = wstrip + pos;
		for (x=xmin; x<xmax; x++, stript++/*, wstript++*/)
		{
			dx = x - d_mx[tid];
			dy = y - d_my[tid];
			if ((r2=dx*dx+dy*dy) < rextlim2)
			{
				if (r2> rintlim2)
				{
					dx += offsetx;
					dy += offsety;
					locarea = 0.0;
					for (sy=APER_OVERSAMP; sy--; dy+=scaley)
					{
						dx1 = dx;
						dy2 = dy*dy;
						for (sx=APER_OVERSAMP; sx--; dx1+=scalex)
							if (dx1*dx1+dy2<raper2)
								locarea += scale2;
					}
				}
				else
					locarea = 1.0;
				area += locarea;
				/*------ Here begin tests for pixel and/or weight overflows. Things are a */
				/*------ bit intricated to have it running as fast as possible in the most */
				/*------ common cases */
				pix=*stript;

				int t = d_labelArray[y*width+x];
				//if ((pix=*stript)<=-BIG || (wfield && (var=*wstript)>=wthresh))
				if(t >=0 && t != tid)
				{
					if (corrflag
							&& (x2=(int)(2*d_mx[tid]+0.49999-x))>=0 && x2<width
							&& (y2=(int)(2*d_my[tid]+0.49999-y))>=fymin && y2<fymax
							&& ((t = d_labelArray[y2*width+x2]) == tid || t == -1))
					{
						pix=*(strip + (pos = (y2%height)*width + x2));
						//if (wfield)
						//{
						//	var = *(wstrip + pos);
						//	if (var>=wthresh)
						//		pix = var = 0.0;
						//}
					}
					else
					{
						pix = 0.0;
						//if (wfield)
						//	var = 0.0;
					}
				}
				if (pflag)
				{
					pix=exp(pix/ngamma);
					sigtv += var*locarea*pix*pix;
				}
				else
					sigtv += var*locarea;
				tv += locarea*pix;

				//if (gainflag && pix>0.0 && gain>0.0)
				//	sigtv += pix/gain*var/backnoise2;
			}
		}
	}

	if (pflag)
	{
		tv = ngamma*(tv-area*exp(d_dbkg[tid]/ngamma));
		sigtv /= ngamma*ngamma;
	}
	else
	{
		tv -= area*d_dbkg[tid];
		if (/*!gainflag && */gain > 0.0 && tv>0.0)
			sigtv += tv/gain;
	}

	if (i<d_prefs->flux_apersize)
		d_flux_aper[i][tid] = tv;
	if (i<d_prefs->fluxerr_apersize)
		d_fluxerr_aper[i][tid] = sqrt(sigtv);
	if (i<d_prefs->mag_apersize)
		d_mag_aper[i][tid] = tv>0.0? -2.5*log10(tv) + d_prefs->mag_zeropoint : 99.0;
	if (i<d_prefs->magerr_apersize)
		d_magerr_aper[i][tid] = tv>0.0? 1.086*sqrt(sigtv)/tv:99.0;

	return;
}

/***************************** computeautoflux********************************/
/*
Compute the total flux within an automatic elliptical aperture.
 */
__device__ void  computeautoflux(

		float		*d_pixelArray,
		int			*d_labelArray, //allocation map
		unsigned int *d_flag,
		float		*d_a,
		float		*d_b,
		float		*d_cxx,
		float		*d_cyy,
		float		*d_cxy,
		float		*d_dbkg,
		float		*d_kronfactor,
		float		*d_flux_auto,
		float		*d_fluxerr_auto,
		float		*d_mag_auto,
		float		*d_magerr_auto,
		double		*d_mx,
		double		*d_my,
		prefstruct	*d_prefs,
		int 		width,
		int			height,
		int			fymin,
		int			fymax,
		double		ngamma,
		double		gain,
		float		backsig,
		float		flagobj2_mag_auto,
		float		flagobj2_magerr_auto,
		int			tid)
{

	double		sigtv, tv, r1, v1,var,backnoise2;
	float		dx,dy, cx2,cy2,cxy, r2,klim2, dxlim, dylim;
	int			area,areab, x,y, x2,y2, xmin,xmax,ymin,ymax,
	pflag, corrflag, pos;
	PIXTYPE		*stript, pix;

	/* Let's initialize some variables */
	//if (!dfield)
	//	dfield = field;
	//if (dwfield)
	//	dwthresh = dwfield->weight_thresh;
	//wstrip = dwstrip = NULL;
	//if (wfield)
	//	wthresh = wfield->weight_thresh;
	//wstript = dwstript = NULL;
	//w = field->width;
	//h = field->stripheight;
	//fymin = field->ymin;
	//fymax = field->ymax;
	//ngamma = field->ngamma;
	//bkg = (double)obj->dbkg;
	//mx = obj->mx;
	//my = obj->my;
	var = backnoise2 = backsig*backsig;
	//gain = field->gain;
	pflag = (d_prefs->detect_type==d_prefs->PHOTO)? 1:0;
	corrflag = (d_prefs->mask_type==d_prefs->MASK_CORRECT);
	//gainflag = wfield && prefs.weightgain_flag;

	/* First step: find the extent of the ellipse (the kron factor r1) */
	/* Clip boundaries in x and y */
	/* We first check that the search ellipse is large enough... */
	if (KRON_NSIG*sqrt(d_a[tid]*d_b[tid])>d_prefs->autoaper[0]/2.0)
	{
		cx2 = d_cxx[tid];
		cy2 = d_cyy[tid];
		cxy = d_cxy[tid];
		dxlim = cx2 - cxy*cxy/(4.0*cy2);
		dxlim = dxlim>0.0 ? KRON_NSIG/sqrt(dxlim) : 0.0;
		dylim = cy2 - cxy*cxy/(4.0*cx2);
		dylim = dylim > 0.0 ? KRON_NSIG/sqrt(dylim) : 0.0;
		klim2 = KRON_NSIG*KRON_NSIG;
	}
	else
		/*-- ...if not, use the circular aperture provided by the user */
	{
		cx2 = cy2 = 1.0;
		cxy = 0.0;
		dxlim = dylim = d_prefs->autoaper[0]/2.0;
		klim2 =  dxlim*dxlim;
	}

	if ((xmin = RINT(d_mx[tid]-dxlim)) < 0)
	{
		xmin = 0;
		d_flag[tid] |= OBJ_APERT_PB;
	}
	if ((xmax = RINT(d_mx[tid]+dxlim)+1) > width)
	{
		xmax = width;
		d_flag[tid] |= OBJ_APERT_PB;
	}
	if ((ymin = RINT(d_my[tid]-dylim)) < fymin)
	{
		ymin = fymin;
		d_flag[tid] |= OBJ_APERT_PB;
	}
	if ((ymax = RINT(d_my[tid]+dylim)+1) > fymax)
	{
		ymax = fymax;
		d_flag[tid] |= OBJ_APERT_PB;
	}

	v1 = r1 = 0.0;
	area = areab = 0;

	//if (dwfield)
	//	dwstrip = dwfield->strip;
	int t;

	for (y=ymin; y<ymax; y++)
	{
		stript = d_pixelArray + (pos = xmin + (y%height)*width);
		//if (dwfield)
		//	dwstript = dwstrip + pos;
		for (x=xmin; x<xmax; x++, stript++)
		{
			dx = x - d_mx[tid];
			dy = y - d_my[tid];
			if ((r2=cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) <= klim2)
			{
				pix=*stript;
				t = d_labelArray[y*width+x];

				if((t == tid || t == -1) & pix>-BIG)
				//if ((pix=*dstript)>-BIG && (!dwfield || (dwfield&&*dwstript<dwthresh)))
				{
					area++;
					r1 += sqrt(r2)*pix;
					v1 += pix;
				}
				else
					areab++;
			}
		}
	}

	area += areab;

	if (area)
	{
		/*-- Go further only if some pixels are available !! */
		if (r1>0.0 && v1>0.0)
		{
			d_kronfactor[tid] = d_prefs->autoparam[0]*r1/v1;
			if (d_kronfactor[tid] < d_prefs->autoparam[1])
				d_kronfactor[tid] = d_prefs->autoparam[1];
		}
		else
			d_kronfactor[tid] = d_prefs->autoparam[1];

		/*-- Flag if the Kron photometry can be strongly affected by neighhours */
		if ((float)areab/area > CROWD_THRESHOLD)
			d_flag[tid] |= OBJ_CROWDED;

		/*-- Second step: integrate within the ellipse */
		/*-- Clip boundaries in x and y (bis) */
		/*-- We first check that the derived ellipse is large enough... */
		if (d_kronfactor[tid]*sqrt(d_a[tid]*d_b[tid])>d_prefs->autoaper[1]/2.0)
		{
			cx2 = d_cxx[tid];
			cy2 = d_cyy[tid];
			cxy = d_cxy[tid];
			dxlim = cx2 - cxy*cxy/(4.0*cy2);
			dxlim = dxlim>0.0 ? d_kronfactor[tid]/sqrt(dxlim) : 0.0;
			dylim = cy2 - cxy*cxy/(4.0*cx2);
			dylim = dylim > 0.0 ? d_kronfactor[tid]/sqrt(dylim) : 0.0;
			klim2 = d_kronfactor[tid]*d_kronfactor[tid];
		}
		else
			/*---- ...if not, use the circular aperture provided by the user */
		{
			cx2 = cy2 = 1.0;
			cxy = 0.0;
			dxlim = dylim = d_prefs->autoaper[1]/2.0;
			klim2 =  dxlim*dxlim;
			d_kronfactor[tid] = 0.0;
		}

		if ((xmin = RINT(d_mx[tid]-dxlim)) < 0)
		{
			xmin = 0;
			d_flag[tid] |= OBJ_APERT_PB;
		}
		if ((xmax = RINT(d_mx[tid]+dxlim)+1) > width)
		{
			xmax = width;
			d_flag[tid] |= OBJ_APERT_PB;
		}
		if ((ymin = RINT(d_my[tid]-dylim)) < fymin)
		{
			ymin = fymin;
			d_flag[tid] |= OBJ_APERT_PB;
		}
		if ((ymax = RINT(d_my[tid]+dylim)+1) > fymax)
		{
			ymax = fymax;
			d_flag[tid] |= OBJ_APERT_PB;
		}

		area = areab = 0;
		tv = sigtv = 0.0;

		//strip = field->strip;

		//if (wfield)
		//	wstrip = wfield->strip;

		for (y=ymin; y<ymax; y++)
		{
			stript = d_pixelArray + (pos = xmin + (y%height)*width);
			//if (wfield)
			//	wstript = wstrip + pos;
			for (x=xmin; x<xmax; x++, stript++)
			{
				dx = x - d_mx[tid];
				dy = y - d_my[tid];
				if ((cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) <= klim2)
				{
					area++;
					/*-------- Here begin tests for pixel and/or weight overflows. Things are a */
					/*-------- bit intricated to have it running as fast as possible in the most */
					/*-------- common cases */
					pix=*stript;

					t = d_labelArray[y*width+x];
					if(t >= 0 && t != tid)
					//if ((pix=*stript)<=-BIG || (wfield && (var=*wstript)>=wthresh))
					{
						areab++;
						if (corrflag
								&& (x2=(int)(2*d_mx[tid]+0.49999-x))>=0 && x2<width
								&& (y2=(int)(2*d_my[tid]+0.49999-y))>=fymin && y2<fymax
								&& ((t = d_labelArray[y2*width+x2]) == -1
										|| t == tid))
						{
							pix=*(d_pixelArray + (pos = (y2%height)*width + x2));
							//if (wfield)
							//{
							//	var = *(wstrip + pos);
							//	if (var>=wthresh)
							//		pix = var = 0.0;
							//}
						}
						else
						{
							pix = 0.0;
							//if (wfield)
							//	var = 0.0;
						}
					}
					if (pflag)
					{
						pix = exp(pix/ngamma);
						sigtv += var*pix*pix;
					}
					else
						sigtv += var;
					tv += pix;
					//if (gainflag && pix>0.0 && gain>0.0)
					//	sigtv += pix/gain*var/backnoise2;
				}
			}
		}

		/*-- Flag if the Kron photometry can be strongly affected by neighhours */
		if ((float)areab > CROWD_THRESHOLD*area)
			d_flag[tid] |= OBJ_CROWDED;

		if (pflag)
		{
			tv = ngamma*(tv-area*exp(((double)d_dbkg[tid])/ngamma));
			sigtv /= ngamma*ngamma;
		}
		else
		{
			tv -= area*(double)d_dbkg[tid];
			if (/*!gainflag &&*/ gain > 0.0 && tv>0.0)
				sigtv += tv/gain;
		}
	}
	else
		/*-- No available pixels: set the flux to zero */
		tv = sigtv = 0.0;

	d_flux_auto[tid] = tv;
	d_fluxerr_auto[tid] = sqrt(sigtv);

	if (*((char *)&flagobj2_mag_auto))
		d_mag_auto[tid] = d_flux_auto[tid]>0.0?
				-2.5*log10(d_flux_auto[tid]) + d_prefs->mag_zeropoint :99.0;
	if (*((char *)&flagobj2_magerr_auto))
		d_magerr_auto[tid] = d_flux_auto[tid]>0.0?
				1.086*d_fluxerr_auto[tid]/d_flux_auto[tid] :99.0;
	if (tv<=0.0)
		d_kronfactor[tid] = 0.0;

	return;
}

/******************************** neurresp **********************************/

/*
Sigmoid function for a neural network.
*/
__device__ double ff(double x)

  {
  return 1.0 / (1.0 + exp(-x));
  }

/*
Neural network response to an input vector.
 */
__device__ void	neurresp(brainstruct *d_brain, double *input, double *output)

{
	brainstruct t_brain = *d_brain;

	int	i, j, l, lastlay = t_brain.layersnb-1;
	double	neursum;

	for (i=0; i<t_brain.nn[0]; i++)
		t_brain.n[0][i] = input[i]*t_brain.inscale[i] + t_brain.inbias[i];

	for (l=0; l<lastlay; l++)
		for (j=0; j<t_brain.nn[l+1]; j++)
		{
			neursum = t_brain.b[l][j];
			for (i=0; i<t_brain.nn[l]; i++)
				neursum += t_brain.w[l][i][j] * t_brain.n[l][i];
			t_brain.n[l+1][j] = ff(neursum);
		}
	for (i=0; i<t_brain.nn[lastlay]; i++)
		output[i] = (t_brain.n[lastlay][i]-t_brain.outbias[i])
		/ t_brain.outscale[i];

	return;
}

/******************************* computemags *********************************/
/*
Compute magnitude parameters.
 */
__device__ void  computemags(
		float	*d_mag_iso,
		float	*d_magerr_iso,
		float	*d_flux_iso,
		float	*d_fluxerr_iso,
		float 	flagobj2_mag_iso,
		float 	flagobj2_magerr_iso,
		double	prefs_mag_zeropoint,
		int tid)

{
	/* Mag. isophotal */
	if (*((char *)&flagobj2_mag_iso)) //true
		d_mag_iso[tid] = d_flux_iso[tid]>0.0?
				-2.5*log10(d_flux_iso[tid]) + prefs_mag_zeropoint:99.0;
	if (*((char *)&flagobj2_magerr_iso))		  //true
		d_magerr_iso[tid] = d_flux_iso[tid]>0.0?
				1.086*d_fluxerr_iso[tid]/d_flux_iso[tid]:99.0;

//	/* Mag. isophotal corrected */
//	if (FLAG(obj2.mag_isocor)) //false
//		obj2->mag_isocor = obj2->flux_isocor>0.0?
//				-2.5*log10(obj2->flux_isocor) + prefs.mag_zeropoint
//				:99.0;
//	if (FLAG(obj2.magerr_isocor)) //false
//		obj2->magerr_isocor = obj2->flux_isocor>0.0?
//				1.086*obj2->fluxerr_isocor/obj2->flux_isocor
//				:99.0;
//
//	/* Choose the ``best'' flux according to the local crowding */
//
//	if (FLAG(obj2.flux_best)) //false
//	{
//		if (obj->flag&OBJ_CROWDED)
//		{
//			obj2->flux_best = obj2->flux_isocor;
//			obj2->fluxerr_best = obj2->fluxerr_isocor;
//		}
//		else
//		{
//			obj2->flux_best = obj2->flux_auto;
//			obj2->fluxerr_best = obj2->fluxerr_auto;
//		}
//	}
//
//	/* Mag. Best */
//	if (FLAG(obj2.mag_best)) //false
//		obj2->mag_best = obj2->flux_best>0.0?
//				-2.5*log10(obj2->flux_best) + prefs.mag_zeropoint
//				:99.0;
//	if (FLAG(obj2.magerr_best)) //false
//		obj2->magerr_best = obj2->flux_best>0.0?
//				1.086*obj2->fluxerr_best/obj2->flux_best
//				:99.0;
//
//	/* Mag. SOM-fit */
//	if (FLAG(obj2.mag_somfit))	//false
//		obj2->mag_somfit = obj2->flux_somfit>0.0?
//				-2.5*log10(obj2->flux_somfit) + prefs.mag_zeropoint
//				:99.0;
//	if (FLAG(obj2.magerr_somfit))	//false
//		obj2->magerr_somfit = obj2->flux_somfit>0.0?
//				1.086*obj2->fluxerr_somfit/obj2->flux_somfit
//				:99.0;
//
//	/* Mag. models */
//	if (FLAG(obj2.mag_prof))	//false
//		obj2->mag_prof = obj2->flux_prof>0.0?
//				-2.5*log10(obj2->flux_prof) + prefs.mag_zeropoint
//				:99.0;
//	if (FLAG(obj2.magerr_prof))		//false
//		obj2->magerr_prof = obj2->flux_prof>0.0?
//				1.086*obj2->fluxerr_prof/obj2->flux_prof
//				:99.0;
//
//	if (FLAG(obj2.prof_spheroid_mag))		//false
//		obj2->prof_spheroid_mag = obj2->prof_spheroid_flux>0.0?
//				-2.5*log10(obj2->prof_spheroid_flux)
//	+ prefs.mag_zeropoint
//	:99.0;
//	if (FLAG(obj2.prof_spheroid_magerr))	//false
//		obj2->prof_spheroid_magerr = obj2->prof_spheroid_flux>0.0?
//				1.086*obj2->prof_spheroid_fluxerr
//				/ obj2->prof_spheroid_flux
//				:99.0;
//
//	if (FLAG(obj2.prof_disk_mag))		//false
//		obj2->prof_disk_mag = obj2->prof_disk_flux>0.0?
//				-2.5*log10(obj2->prof_disk_flux)
//	+ prefs.mag_zeropoint
//	:99.0;
//	if (FLAG(obj2.prof_disk_magerr))	//false
//		obj2->prof_disk_magerr = obj2->prof_disk_flux>0.0?
//				1.086*obj2->prof_disk_fluxerr
//				/ obj2->prof_disk_flux
//				:99.0;
//
//	if (FLAG(obj2.prof_bar_mag))		//false
//		obj2->prof_bar_mag = obj2->prof_bar_flux>0.0?
//				-2.5*log10(obj2->prof_bar_flux)
//	+ prefs.mag_zeropoint
//	:99.0;
//	if (FLAG(obj2.prof_bar_magerr))		//false
//		obj2->prof_bar_magerr = obj2->prof_bar_flux>0.0?
//				1.086*obj2->prof_bar_fluxerr
//				/obj2->prof_bar_flux
//				:99.0;
//
//	if (FLAG(obj2.prof_arms_mag))		//false
//		obj2->prof_arms_mag = obj2->prof_arms_flux>0.0?
//				-2.5*log10(obj2->prof_arms_flux)
//	+ prefs.mag_zeropoint
//	:99.0;
//	if (FLAG(obj2.prof_arms_magerr))	//false
//		obj2->prof_arms_magerr = obj2->prof_arms_flux>0.0?
//				1.086*obj2->prof_arms_fluxerr
//				/obj2->prof_arms_flux
//				:99.0;

//	/* Mag. WINdowed */
//	if (FLAG(obj2.mag_win))		//false
//		obj2->mag_win = obj2->flux_win>0.0?
//				-2.5*log10(obj2->flux_win) + prefs.mag_zeropoint
//				:99.0;
//	if (FLAG(obj2.magerr_win))	//false
//		obj2->magerr_win = obj2->flux_win>0.0?
//				1.086*obj2->fluxerr_win/obj2->flux_win
//				:99.0;
//	/* Mag. GALFIT */
//	if (FLAG(obj2.mag_galfit))	//false
//		obj2->mag_galfit = obj2->flux_galfit>0.0?
//				-2.5*log10(obj2->flux_galfit) + prefs.mag_zeropoint
//				:99.0;
//	if (FLAG(obj2.magerr_galfit)) //false
//		obj2->magerr_galfit = obj2->flux_galfit>0.0?
//				1.086*obj2->fluxerr_galfit/obj2->flux_galfit
//				:99.0;

//	/* SB units */
//	if (FLAG(obj2.maxmu))	//false
//		outobj2.maxmu = obj->peak > 0.0 ?
//				-2.5*log10((obj->peak)
//						/ (field->pixscale * field->pixscale)) + prefs.mag_zeropoint
//						: 99.0;
//
//	if (FLAG(obj2.threshmu))	//false
//		obj2->threshmu = obj->thresh > 0.0 ?
//				-2.5*log10((obj->thresh)
//						/ (field->pixscale * field->pixscale)) + prefs.mag_zeropoint:99.0;

	return;
}

/*--------------------------------------------------------------------------*/

//__device__ int  npcode = 26;
//__device__ char pcodes[26][4] =
//      {"AZP", "TAN", "SIN", "STG", "ARC", "ZPN", "ZEA", "AIR", "CYP", "CAR",
//       "MER", "CEA", "COP", "COD", "COE", "COO", "BON", "PCO", "GLS", "PAR",
//       "AIT", "MOL", "CSC", "QSC", "TSC", "TNX"};
//
//__device__ int wcsset (const int naxis, const char ctype[][9], wcsprm *wcs)
//{
//   int j, k, *ndx;
//   char requir[9];
//
//   strcpy(wcs->pcode, "");
//   strcpy(requir, "");
//   wcs->lng = -1;
//   ndx = &wcs->lng;	/* to satisfy gcc -Wall */
//   wcs->lat = -1;
//   wcs->cubeface = -1;
//
//   for (j = 0; j < naxis; j++) {
//      if (ctype[j][4] != '-') {
//         if (strcmp(ctype[j], "CUBEFACE") == 0) {
//            if (wcs->cubeface == -1) {
//               wcs->cubeface = j;
//            } else {
//               /* Multiple CUBEFACE axes! */
//               return 1;
//            }
//         }
//         continue;
//      }
//
//      /* Got an axis qualifier, is it a recognized WCS projection? */
//      for (k = 0; k < npcode; k++) {
//         if (strncmp(&ctype[j][5], pcodes[k], 3) == 0) break;
//      }
//
//      if (k == npcode) {
//         /* Allow NCP to pass (will be converted to SIN later). */
//         if (strncmp(&ctype[j][5], "NCP", 3)) continue;
//      }
//
//      /* Parse the celestial axis type. */
//      if (strcmp(wcs->pcode, "") == 0) {
//         sprintf(wcs->pcode, "%.3s", &ctype[j][5]);
//
//         if (strncmp(ctype[j], "RA--", 4) == 0) {
//            wcs->lng = j;
//            strcpy(wcs->lngtyp, "RA");
//            strcpy(wcs->lattyp, "DEC");
//            ndx = &wcs->lat;
//            sprintf(requir, "DEC--%s", wcs->pcode);
//         } else if (strncmp(ctype[j], "DEC-", 4) == 0) {
//            wcs->lat = j;
//            strcpy(wcs->lngtyp, "RA");
//            strcpy(wcs->lattyp, "DEC");
//            ndx = &wcs->lng;
//            sprintf(requir, "RA---%s", wcs->pcode);
//         } else if (strncmp(&ctype[j][1], "LON", 3) == 0) {
//            wcs->lng = j;
//            sprintf(wcs->lngtyp, "%cLON", ctype[j][0]);
//            sprintf(wcs->lattyp, "%cLAT", ctype[j][0]);
//            ndx = &wcs->lat;
//            sprintf(requir, "%s-%s", wcs->lattyp, wcs->pcode);
//         } else if (strncmp(&ctype[j][1], "LAT", 3) == 0) {
//            wcs->lat = j;
//            sprintf(wcs->lngtyp, "%cLON", ctype[j][0]);
//            sprintf(wcs->lattyp, "%cLAT", ctype[j][0]);
//            ndx = &wcs->lng;
//            sprintf(requir, "%s-%s", wcs->lngtyp, wcs->pcode);
//         } else {
//            /* Unrecognized celestial type. */
//            return 1;
//         }
//      } else {
//         if (strncmp(ctype[j], requir, 8) != 0) {
//            /* Inconsistent projection types. */
//            return 1;
//         }
//
//         *ndx = j;
//         strcpy(requir, "");
//       }
//     }
//
//   if (strcmp(requir, "")) {
//      /* Unmatched celestial axis. */
//      return 1;
//   }
//
//   if (strcmp(wcs->pcode, "")) {
//      wcs->flag = WCSSET;
//   } else {
//      /* Signal for no celestial axis pair. */
//      wcs->flag = 999;
//   }
//
//   return 0;
//}
//
//__device__ int matinv(int n, double *mat, double *inv)
//{
//   register int i, ij, ik, j, k, kj, pj;
//   int    itemp, mem, *mxl, *lxm, pivot;
//   double colmax, *lu, *rowmax, dtemp;
//
//   /* Allocate memory for internal arrays. */
//   mem = n * sizeof(int);
//   if ((mxl = (int*)malloc(mem)) == (int*)0) return 1;
//   if ((lxm = (int*)malloc(mem)) == (int*)0) {
//      free(mxl);
//      return 1;
//   }
//
//   mem = n * sizeof(double);
//   if ((rowmax = (double*)malloc(mem)) == (double*)0) {
//      free(mxl);
//      free(lxm);
//      return 1;
//   }
//
//   mem *= n;
//   if ((lu = (double*)malloc(mem)) == (double*)0) {
//      free(mxl);
//      free(lxm);
//      free(rowmax);
//      return 1;
//   }
//
//   /* Initialize arrays. */
//   for (i = 0, ij = 0; i < n; i++) {
//      /* Vector which records row interchanges. */
//      mxl[i] = i;
//
//      rowmax[i] = 0.0;
//
//      for (j = 0; j < n; j++, ij++) {
//         dtemp = fabs(mat[ij]);
//         if (dtemp > rowmax[i]) rowmax[i] = dtemp;
//
//         lu[ij] = mat[ij];
//      }
//
//      /* A row of zeroes indicates a singular matrix. */
//      if (rowmax[i] == 0.0) {
//         free(mxl);
//         free(lxm);
//         free(rowmax);
//         free(lu);
//         return 2;
//      }
//   }
//
//
//   /* Form the LU triangular factorization using scaled partial pivoting. */
//   for (k = 0; k < n; k++) {
//      /* Decide whether to pivot. */
//      colmax = fabs(lu[k*n+k]) / rowmax[k];
//      pivot = k;
//
//      for (i = k+1; i < n; i++) {
//         ik = i*n + k;
//         dtemp = fabs(lu[ik]) / rowmax[i];
//         if (dtemp > colmax) {
//            colmax = dtemp;
//            pivot = i;
//         }
//      }
//
//      if (pivot > k) {
//         /* We must pivot, interchange the rows of the design matrix. */
//         for (j = 0, pj = pivot*n, kj = k*n; j < n; j++, pj++, kj++) {
//            dtemp = lu[pj];
//            lu[pj] = lu[kj];
//            lu[kj] = dtemp;
//         }
//
//         /* Amend the vector of row maxima. */
//         dtemp = rowmax[pivot];
//         rowmax[pivot] = rowmax[k];
//         rowmax[k] = dtemp;
//
//         /* Record the interchange for later use. */
//         itemp = mxl[pivot];
//         mxl[pivot] = mxl[k];
//         mxl[k] = itemp;
//      }
//
//      /* Gaussian elimination. */
//      for (i = k+1; i < n; i++) {
//         ik = i*n + k;
//
//         /* Nothing to do if lu[ik] is zero. */
//         if (lu[ik] != 0.0) {
//            /* Save the scaling factor. */
//            lu[ik] /= lu[k*n+k];
//
//            /* Subtract rows. */
//            for (j = k+1; j < n; j++) {
//               lu[i*n+j] -= lu[ik]*lu[k*n+j];
//            }
//         }
//      }
//   }
//
//
//   /* mxl[i] records which row of mat corresponds to row i of lu.  */
//   /* lxm[i] records which row of lu  corresponds to row i of mat. */
//   for (i = 0; i < n; i++) {
//      lxm[mxl[i]] = i;
//   }
//
//
//   /* Determine the inverse matrix. */
//   for (i = 0, ij = 0; i < n; i++) {
//      for (j = 0; j < n; j++, ij++) {
//         inv[ij] = 0.0;
//      }
//   }
//
//   for (k = 0; k < n; k++) {
//      inv[lxm[k]*n+k] = 1.0;
//
//      /* Forward substitution. */
//      for (i = lxm[k]+1; i < n; i++) {
//         for (j = lxm[k]; j < i; j++) {
//            inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
//         }
//      }
//
//      /* Backward substitution. */
//      for (i = n-1; i >= 0; i--) {
//         for (j = i+1; j < n; j++) {
//            inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
//         }
//         inv[i*n+k] /= lu[i*n+i];
//      }
//   }
//
//   free(mxl);
//   free(lxm);
//   free(rowmax);
//   free(lu);
//
//   return 0;
//}
//
//__device__ int linset(linprm *lin)
//{
//   int i, ij, j, mem, n;
//
//   n = lin->naxis;
//
//   /* Allocate memory for internal arrays. */
//   mem = n * n * sizeof(double);
//   lin->piximg = (double*)malloc(mem);
//   if (lin->piximg == (double*)0) return 1;
//
//   lin->imgpix = (double*)malloc(mem);
//   if (lin->imgpix == (double*)0) {
//      free(lin->piximg);
//      return 1;
//   }
//
//   /* Compute the pixel-to-image transformation matrix. */
//   for (i = 0, ij = 0; i < n; i++) {
//      for (j = 0; j < n; j++, ij++) {
//         lin->piximg[ij] = lin->cdelt[i] * lin->pc[ij];
//      }
//   }
//
//   /* Compute the image-to-pixel transformation matrix. */
//   if (matinv(n, lin->piximg, lin->imgpix)) return 2;
//
//   lin->flag = LINSET;
//
//   return 0;
//}
//
//__device__ int linrev(const double pixcrd[], linprm *lin, double imgcrd[])
//{
//   int i, ij, j, n;
//   double temp;
//
//   n = lin->naxis;
//
////   if (lin->flag != LINSET) {
////      if (linset(lin)) return 1;
////   }
//
//   for (i = 0; i < n; i++) {
//      imgcrd[i] = 0.0;
//   }
//
//   for (j = 0; j < n; j++) {
//      temp = pixcrd[j] - lin->crpix[j];
//      for (i = 0, ij = j; i < n; i++, ij+=n) {
//         imgcrd[i] += lin->piximg[ij] * temp;
//      }
//   }
//
//   return 0;
//}
//
//__device__ double wcs_cosd(const double angle)
//{
//   double resid;
//
//   const double d2r = PI / 180.0;
//
//   resid = fabs(fmod(angle,360.0));
//   if (resid == 0.0) {
//      return 1.0;
//   } else if (resid == 90.0) {
//      return 0.0;
//   } else if (resid == 180.0) {
//      return -1.0;
//   } else if (resid == 270.0) {
//      return 0.0;
//   }
//
//   return cos(angle*d2r);
//}
//
//__device__ double wcs_sind(const double angle)
//{
//   double resid;
//
//   const double d2r = PI / 180.0;
//
//   resid = fmod(angle-90.0,360.0);
//   if (resid == 0.0) {
//      return 1.0;
//   } else if (resid == 90.0) {
//      return 0.0;
//   } else if (resid == 180.0) {
//      return -1.0;
//   } else if (resid == 270.0) {
//      return 0.0;
//   }
//
//   return sin(angle*d2r);
//}
//
//__device__ double wcs_atan2d(double y, double x)
//{
//	const double r2d = 180.0 / PI;
//
//   if (y == 0.0) {
//      if (x >= 0.0) {
//         return 0.0;
//      } else if (x < 0.0) {
//         return 180.0;
//      }
//   } else if (x == 0.0) {
//      if (y > 0.0) {
//         return 90.0;
//      } else if (y < 0.0) {
//         return -90.0;
//      }
//   }
//
//   return atan2(y,x)*r2d;
//}
//
//__device__ double wcs_acosd(double v)
//{
//	const double r2d = 180.0 / PI;
//
//   if (v >= 1.0) {
//      if (v-1.0 <  WCSTRIG_TOL) return 0.0;
//   } else if (v == 0.0) {
//      return 90.0;
//   } else if (v <= -1.0) {
//      if (v+1.0 > -WCSTRIG_TOL) return 180.0;
//   }
//
//   return acos(v)*r2d;
//}
//
//__device__ double wcs_asind(double v)
//{
//	const double r2d = 180.0 / PI;
//
//   if (v <= -1.0) {
//      if (v+1.0 > -WCSTRIG_TOL) return -90.0;
//   } else if (v == 0.0) {
//      return 0.0;
//   } else if (v >= 1.0) {
//      if (v-1.0 <  WCSTRIG_TOL) return 90.0;
//   }
//
//   return asin(v)*r2d;
//}
//
//#define wcs_copysign(X, Y) ((Y) < 0.0 ? -fabs(X) : fabs(X))
//
//__device__ int sphrev (double phi, double theta, double *eul, double *lng, double *lat)
//{
//   double cosphi, costhe, dlng, dphi, sinphi, sinthe, x, y, z;
//   const double tol = 1.0e-5;
//
//   costhe = wcs_cosd(theta);
//   sinthe = wcs_sind(theta);
//
//   dphi = phi - eul[2];
//   cosphi = wcs_cosd(dphi);
//   sinphi = wcs_sind(dphi);
//
//   /* Compute the celestial longitude. */
//   x = sinthe*eul[4] - costhe*eul[3]*cosphi;
//   if (fabs(x) < tol) {
//      /* Rearrange formula to reduce roundoff errors. */
//      x = -wcs_cosd(theta+eul[1]) + costhe*eul[3]*(1.0 - cosphi);
//   }
//   y = -costhe*sinphi;
//   if (x != 0.0 || y != 0.0) {
//      dlng = wcs_atan2d(y, x);
//   } else {
//      /* Change of origin of longitude. */
//      dlng = dphi + 180.0;
//   }
//   *lng = eul[0] + dlng;
//
//   /* Normalize the celestial longitude. */
//   if (eul[0] >= 0.0) {
//      if (*lng < 0.0) *lng += 360.0;
//   } else {
//      if (*lng > 0.0) *lng -= 360.0;
//   }
//
//   if (*lng > 360.0) {
//      *lng -= 360.0;
//   } else if (*lng < -360.0) {
//      *lng += 360.0;
//   }
//
//   /* Compute the celestial latitude. */
//   if (fmod(dphi,180.0) == 0.0) {
//      *lat = theta + cosphi*eul[1];
//      if (*lat >  90.0) *lat =  180.0 - *lat;
//      if (*lat < -90.0) *lat = -180.0 - *lat;
//   } else {
//      z = sinthe*eul[3] + costhe*eul[4]*cosphi;
//      if (fabs(z) > 0.99) {
//         /* Use an alternative formula for greater numerical accuracy. */
//         *lat = wcs_copysign(wcs_acosd(sqrt(x*x+y*y)), z);
//      } else {
//         *lat = wcs_asind(z);
//      }
//   }
//
//   return 0;
//}
//
//__device__ int celrev(const char pcode[4],
//		const double x,
//		const double y,
//		prjprm *prj,
//		double *phi,
//		double *theta,
//		celprm *cel,
//		double *lng,
//		double *lat)
//{
//   int    err;
//
////   if (cel->flag != CELSET) {
////      if(celset(pcode, cel, prj)) return 1;
////   }
//
//   /* Apply reverse projection. */
//   if ((err = cel->prjrev(x, y, prj, phi, theta))) {
//      return err == 1 ? 2 : 3;
//   }
//   if (fabs(*phi)>180.0 || fabs(*theta)>90.0)
//     return 2;
//
//   /* Compute native coordinates. */
//   sphrev(*phi, *theta, cel->euler, lng, lat);
//
//   return 0;
//}
//
///*--------------------------------------------------------------------------*/
//__device__ int wcsrev(const char ctype[][9],
//		wcsprm *wcs,
//		const double pixcrd[],
//		linprm *lin,
//		double imgcrd[],
//		prjprm *prj,
//		double *phi,
//		double *theta,
//		const double crval[],
//		celprm *cel,
//		double world[])
//{
//   int    err, face, j;
//   double offset;
//
////   /* Initialize if required. */
////   if (wcs->flag != WCSSET) {
////      if (wcsset(lin->naxis, ctype, wcs)) return 1;
////   }
//
//   /* Apply reverse linear transformation. */
//   if (linrev(pixcrd, lin, imgcrd)) {
//      return 4;
//   }
//
//   /* Convert to world coordinates. */
//   for (j = 0; j < lin->naxis; j++) {
//      if (j == wcs->lng) continue;
//      if (j == wcs->lat) continue;
//      world[j] = imgcrd[j] + crval[j];
//   }
//
//
//   if (wcs->flag != 999) {
//      /* Do we have a CUBEFACE axis? */
//      if (wcs->cubeface != -1) {
//         face = (int)(imgcrd[wcs->cubeface] + 0.5);
//         if (fabs(imgcrd[wcs->cubeface]-face) > 1e-10) {
//            return 3;
//         }
//
//         /* Separation between faces. */
//         if (prj->r0 == 0.0) {
//            offset = 90.0;
//         } else {
//            offset = prj->r0*PI/2.0;
//         }
//
//         /* Lay out faces in a plane. */
//         switch (face) {
//         case 0:
//            imgcrd[wcs->lat] += offset;
//            break;
//         case 1:
//            break;
//         case 2:
//            imgcrd[wcs->lng] += offset;
//            break;
//         case 3:
//            imgcrd[wcs->lng] += offset*2;
//            break;
//         case 4:
//            imgcrd[wcs->lng] += offset*3;
//            break;
//         case 5:
//            imgcrd[wcs->lat] -= offset;
//            break;
//         default:
//            return 3;
//         }
//      }
//
//      /* Compute celestial coordinates. */
//      if (strcmp(wcs->pcode, "NCP") == 0) {
//         /* Convert NCP to SIN. */
//         if (cel->ref[2] == 0.0) {
//            return 2;
//         }
//
//         strcpy(wcs->pcode, "SIN");
//         prj->p[1] = 0.0;
//         prj->p[2] = wcs_cosd(cel->ref[2])/wcs_sind(cel->ref[2]);
//         prj->flag = 0;
//      }
//
//      if ((err = celrev(wcs->pcode, imgcrd[wcs->lng], imgcrd[wcs->lat], prj,
//                   phi, theta, cel, &world[wcs->lng], &world[wcs->lat]))) {
//         return err;
//      }
//   }
//
//   return 0;
//}
//
///*--------------------------------------------------------------------------*/
//__device__ int	celsys_to_eq(wcsstruct *wcs, double *wcspos)
//
//{
//	double	*mat,
//	a2,d2,sd2,cd2cp,sd,x,y;
//	int		lng, lat;
//
//	mat = wcs->celsysmat;
//	a2 = wcspos[lng = wcs->wcsprm->lng]*DEG - mat[1];
//	d2 = wcspos[lat = wcs->wcsprm->lat]*DEG;
//	/* A bit of spherical trigonometry... */
//	/* Compute the latitude... */
//	sd2 = sin(d2);
//	cd2cp = cos(d2)*mat[2];
//	sd = sd2*mat[3]-cd2cp*cos(a2);
//	/* ...and the longitude */
//	y = cd2cp*sin(a2);
//	x = sd2 - sd*mat[3];
//	wcspos[lng] = fmod((atan2(y,x) + mat[0])/DEG+360.0, 360.0);
//	wcspos[lat] = asin(sd)/DEG;
//
//	return RETURN_OK;
//}
//
///******* raw_to_wcs ***********************************************************
//PROTO	int raw_to_wcs(wcsstruct *, double *, double *)
//PURPOSE	Convert raw (pixel) coordinates to WCS (World Coordinate System).
//INPUT	WCS structure,
//	Pointer to the array of input coordinates,
//	Pointer to the array of output coordinates.
//OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
//NOTES	-.
//AUTHOR	E. Bertin (IAP)
//VERSION	08/02/2007
// ***/
//__device__ int	raw_to_wcs(wcsstruct *d_wcs, double *pixpos, double *wcspos)
//
//{
//	double	imgcrd[NAXIS], phi, theta;
//	int		i;
//
//	if (wcsrev((const char(*)[9])d_wcs->ctype, d_wcs->wcsprm, pixpos,
//			d_wcs->lin,imgcrd, d_wcs->prj, &phi, &theta, d_wcs->crval, d_wcs->cel, wcspos))
//	{
//		for (i=0; i<d_wcs->naxis; i++)
//			wcspos[i] = WCS_NOCOORD;
//		return RETURN_ERROR;
//	}
//
//	/* If needed, convert from a different coordinate system to equatorial */
//	if (d_wcs->celsysconvflag)
//		celsys_to_eq(d_wcs, wcspos);
//
//	return RETURN_OK;
//}
//
//__device__ void	precess(double yearin, double alphain, double deltain,
//		double yearout, double *alphaout, double *deltaout)
//
//{
//	double	dzeta,theta,z, t1,t1t1, t2,t2t2,t2t2t2,
//	cddsadz, cddcadz, cdd, sdd, adz, cdin,sdin,ct,st,caindz;
//
//	alphain *= DEG;
//	deltain *= DEG;
//
//	t1 = (yearin - 2000.0)/1000.0;
//	t2 = (yearout - yearin)/1000.0;
//	t1t1 = t1*t1;
//	t2t2t2 = (t2t2 = t2*t2)*t2;
//	theta = (97171.735e-06 - 413.691e-06*t1 - 1.052e-06 * t1t1) * t2
//			+ (-206.846e-06 - 1.052e-06*t1) * t2t2 - 202.812e-06 * t2t2t2;
//	dzeta = (111808.609e-06 + 677.071e-06*t1 - 0.674e-06 * t1t1) * t2
//			+ (146.356e-06 - 1.673e-06*t1) * t2t2 + 87.257e-06 * t2t2t2;
//	z = (111808.609e-06 +677.071e-06*t1 - 0.674e-06 * t1t1) * t2
//			+ (530.716e-06 + 0.320e-06*t1) * t2t2 + 88.251e-06 * t2t2t2;
//	cddsadz = (cdin=cos(deltain)) * sin(alphain+dzeta);
//	cddcadz = -(sdin=sin(deltain))*(st=sin(theta))
//			+cdin*(ct=cos(theta))*(caindz=cos(alphain+dzeta));
//	sdd = sdin*ct + cdin*st*caindz;
//	cdd = cos(*deltaout = asin(sdd));
//	adz = asin(cddsadz/cdd);
//	if (cddcadz<0.0)
//		adz = PI - adz;
//	if (adz<0.0)
//		adz += 2.0*PI;
//	adz += z;
//	*alphaout = adz/DEG;
//	*deltaout /= DEG;
//
//	return;
//}

/**************************** computeastrom *********************************/
/*
Compute real WORLD coordinates and dimensions according to FITS info.
 */
__device__ void	computeastrom1(
		double 		*d_posx,
		double		*d_posy,
		double		*d_mxw,
		double		*d_myw,
		double		*d_alphas,
		double		*d_deltas,
		double		*d_alpha2000,
		double		*d_delta2000,
		wcsstruct	*d_wcs,
		double		flagobj2_alpha2000,
		double		flagobj2_mxw,
		int			tid)

{
	double	rawpos[NAXIS], wcspos[NAXIS],
	pixscale2, da,dd;
	int		lng,lat;

	//wcs = field->wcs;
	lng = d_wcs->lng;
	lat = d_wcs->lat;
	pixscale2 = 0.0;	/* To avoid gcc -Wall warnings */

	/* If working with WCS, compute WORLD coordinates and local matrix */
	if (*((char *)&flagobj2_mxw))		//true
	{
		rawpos[0] = d_posx[tid];
		rawpos[1] = d_posy[tid];
		//raw_to_wcs(d_wcs, rawpos, wcspos);
		d_mxw[tid] = wcspos[0];
		d_myw[tid] = wcspos[1];
		if (lng != lat)
		{
			d_alphas[tid] = lng<lat? d_mxw[tid] : d_myw[tid];
			d_deltas[tid] = lng<lat? d_myw[tid] : d_mxw[tid];
			if(*((char *)&flagobj2_alpha2000))	//true
			{
				if (fabs(d_wcs->equinox-2000.0)>0.003)
					;//precess(d_wcs->equinox, wcspos[lng<lat?0:1], wcspos[lng<lat?1:0],
					//		2000.0, &d_alpha2000[tid], &d_delta2000[tid]);
				else
				{
					d_alpha2000[tid] = lng<lat? d_mxw[tid] : d_myw[tid];
					d_delta2000[tid] = lng<lat? d_myw[tid] : d_mxw[tid];
				}

				//////////////////////////////////////////////////////
//				if (FLAG(obj2.dtheta2000))	//false
//				{
//					da = wcs->ap2000 - obj2->alpha2000;
//					dd = (sin(wcs->dp2000*DEG)
//							-sin(obj2->delta2000*DEG)*sin(obj2->deltas*DEG))
//							/(cos(obj2->delta2000*DEG)*cos(obj2->deltas*DEG));
//					dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
//					obj2->dtheta2000 = (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
//				}
//				if (FLAG(obj2.alpha1950))	//false
//				{
//					j2b(wcs->equinox, obj2->alpha2000, obj2->delta2000,
//							&obj2->alpha1950, &obj2->delta1950);
//					if (FLAG(obj2.dtheta1950))
//					{
//						da = wcs->ap1950 - obj2->alpha1950;
//						dd = (sin(wcs->dp1950*DEG)
//								-sin(obj2->delta1950*DEG)*sin(obj2->deltas*DEG))
//								/(cos(obj2->delta1950*DEG)*cos(obj2->deltas*DEG));
//						dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
//						obj2->dtheta1950 = (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
//					}
//				}
			}
		}
	}

//	/* Idem for peak-flux positions */
//	if (FLAG(obj2.peakxw))	//false
//	{
//		rawpos[0] = obj->peakx;
//		rawpos[1] = obj->peaky;
//		raw_to_wcs(wcs, rawpos, wcspos);
//		obj2->peakxw = wcspos[0];
//		obj2->peakyw = wcspos[1];
//		if (lng != lat)
//		{
//			obj2->peakalphas = lng<lat? obj2->peakxw : obj2->peakyw;
//			obj2->peakdeltas = lng<lat? obj2->peakyw : obj2->peakxw;
//			if (FLAG(obj2.peakalpha2000))
//			{
//				if (fabs(wcs->equinox-2000.0)>0.003)
//					precess(wcs->equinox, wcspos[lng<lat?0:1], wcspos[lng<lat?1:0],
//							2000.0, &obj2->peakalpha2000, &obj2->peakdelta2000);
//				else
//				{
//					obj2->peakalpha2000 = lng<lat? obj2->peakxw : obj2->peakyw;
//					obj2->peakdelta2000 = lng<lat? obj2->peakyw : obj2->peakxw;
//				}
//				if (FLAG(obj2.peakalpha1950))
//					j2b(wcs->equinox, obj2->peakalpha2000, obj2->peakdelta2000,
//							&obj2->peakalpha1950, &obj2->peakdelta1950);
//			}
//		}
//	}
//
//	/* Idem for Windowed positions */
//	if (FLAG(obj2.winpos_xw))	//false
//	{
//		rawpos[0] = obj2->winpos_x;
//		rawpos[1] = obj2->winpos_y;
//		raw_to_wcs(wcs, rawpos, wcspos);
//		obj2->winpos_xw = wcspos[0];
//		obj2->winpos_yw = wcspos[1];
//		if (lng != lat)
//		{
//			obj2->winpos_alphas = lng<lat? obj2->winpos_xw : obj2->winpos_yw;
//			obj2->winpos_deltas = lng<lat? obj2->winpos_yw : obj2->winpos_xw;
//			if (FLAG(obj2.winpos_alpha2000))
//			{
//				if (fabs(wcs->equinox-2000.0)>0.003)
//					precess(wcs->equinox, wcspos[0], wcspos[1],
//							2000.0, &obj2->winpos_alpha2000, &obj2->winpos_delta2000);
//				else
//				{
//					obj2->winpos_alpha2000 = lng<lat? obj2->winpos_xw : obj2->winpos_yw;
//					obj2->winpos_delta2000 = lng<lat? obj2->winpos_yw : obj2->winpos_xw;
//				}
//				if (FLAG(obj2.winpos_alpha1950))
//					j2b(wcs->equinox, obj2->winpos_alpha2000, obj2->winpos_delta2000,
//							&obj2->winpos_alpha1950, &obj2->winpos_delta1950);
//			}
//		}
//	}
//
//	/* Idem for Model-fitted positions */
//	if (FLAG(obj2.xw_prof))		//false
//	{
//		rawpos[0] = obj2->x_prof;
//		rawpos[1] = obj2->y_prof;
//		raw_to_wcs(wcs, rawpos, wcspos);
//		obj2->xw_prof = wcspos[0];
//		obj2->yw_prof = wcspos[1];
//		if (lng != lat)
//		{
//			obj2->alphas_prof = lng<lat? obj2->xw_prof : obj2->yw_prof;
//			obj2->deltas_prof = lng<lat? obj2->yw_prof : obj2->xw_prof;
//			if (FLAG(obj2.alpha2000_prof))
//			{
//				if (fabs(wcs->equinox-2000.0)>0.003)
//					precess(wcs->equinox, wcspos[0], wcspos[1],
//							2000.0, &obj2->alpha2000_prof, &obj2->delta2000_prof);
//				else
//				{
//					obj2->alpha2000_prof = lng<lat? obj2->xw_prof : obj2->yw_prof;
//					obj2->delta2000_prof = lng<lat? obj2->yw_prof : obj2->xw_prof;
//				}
//				if (FLAG(obj2.alpha1950_prof))
//					j2b(wcs->equinox, obj2->alpha2000_prof, obj2->delta2000_prof,
//							&obj2->alpha1950_prof, &obj2->delta1950_prof);
//			}
//		}
//	}
//
//	/* Custom coordinate system for the MAMA machine */
//	if (FLAG(obj2.mamaposx))	//false
//	{
//		rawpos[0] = obj2->posx - 0.5;
//		rawpos[1] = obj2->posy - 0.5;
//		raw_to_wcs(wcs, rawpos, wcspos);
//		obj2->mamaposx = wcspos[1]*(MAMA_CORFLEX+1.0);
//		obj2->mamaposy = wcspos[0]*(MAMA_CORFLEX+1.0);
//	}
//
//	if (FLAG(obj2.mx2w)
//			|| FLAG(obj2.win_mx2w)
//			|| FLAG(obj2.poserr_mx2w)
//			|| FLAG(obj2.winposerr_mx2w)
//			|| FLAG(obj2.poserrmx2w_prof)
//			|| FLAG(obj2.prof_flagw)
//			|| ((!prefs.pixel_scale) && (FLAG(obj2.npixw)
//					|| FLAG(obj2.fdnpixw)
//					|| FLAG(obj2.fwhmw))))	//false
//	{
//		rawpos[0] = obj2->posx;
//		rawpos[1] = obj2->posy;
//		pixscale2 = wcs_jacobian(wcs, rawpos, obj2->jacob);
//	}
//
//	/* Express shape parameters in WORLD frame */
//	if (FLAG(obj2.mx2w))		//false
//		astrom_shapeparam(field, obj);
//	if (FLAG(obj2.win_mx2w))	//false
//		astrom_winshapeparam(field, obj);
//	if (FLAG(obj2.prof_flagw))	//false
//		astrom_profshapeparam(field, obj);
//
//	/* Express position error parameters in WORLD frame */
//	if (FLAG(obj2.poserr_mx2w))		//false
//		astrom_errparam(field, obj);
//	if (FLAG(obj2.winposerr_mx2w))	//false
//		astrom_winerrparam(field, obj);
//	if (FLAG(obj2.poserrmx2w_prof))	//false
//		astrom_proferrparam(field, obj);
//
//	if (FLAG(obj2.npixw))	//false
//		obj2->npixw = obj->npix * (prefs.pixel_scale?
//				field->pixscale/3600.0*field->pixscale/3600.0 : pixscale2);
//	if (FLAG(obj2.fdnpixw))	//false
//		obj2->fdnpixw = obj->fdnpix * (prefs.pixel_scale?
//				field->pixscale/3600.0*field->pixscale/3600.0 : pixscale2);
//
//	if (FLAG(obj2.fwhmw))	//false
//		obj2->fwhmw = obj->fwhm * (prefs.pixel_scale?
//				field->pixscale/3600.0 : sqrt(pixscale2));

	return;
}

__global__ void preanalyse_full_kernel(
		unsigned int *d_cuttedObjLevelArray,
		unsigned int *d_cuttedObjIndexArray,
		unsigned int *d_cuttedPixCountArray,
		unsigned int *d_cuttedObjFlagArray,
		float		 *d_cuttedDthreshArray,
		float		 *d_cdPixArray,
		float		 *d_pixelArray,
		int			 *d_labelArray,
		unsigned int **d_pixelIndexArray,
		/* output starts */
		unsigned int *d_xmin,
		unsigned int *d_xmax,
		unsigned int *d_ymin,
		unsigned int *d_ymax,
		unsigned int *d_dnpix,
		double *d_mx,
		double *d_my,
		double *d_mx2,
		double *d_my2,
		double *d_mxy,
		float *d_cxx,
		float *d_cxy,
		float *d_cyy,
		float *d_a,
		float *d_b,
		float *d_theta,
		float *d_abcor,
		float *d_fdpeak,
		float *d_dpeak,
		float *d_fdflux,
		float *d_dflux,
		char  *d_singuflag,
		float *d_amp,
		/* output ends */
		int width,
		int height,
		float thresh,	//900.***
		int plistexist_dthresh,
		int analyse_type,
		unsigned int numobj)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= numobj)
		return;

	//get the level the object belongs to.
	unsigned int level = d_cuttedObjLevelArray[tid];
	//get the index of the object in original pix index array.
	unsigned int index = d_cuttedObjIndexArray[tid];
	//get the number of pixels of the object
	unsigned int pixcount = d_cuttedPixCountArray[tid];

	unsigned int flag = d_cuttedObjFlagArray[tid];

	//get the right pixel index array
	unsigned int *pixIndexArray = d_pixelIndexArray[level];

	//start to compute the attributes
	/*-----  initialize stacks and bounds */
	float	peak, cpeak, val, cval, minthresh, thresht;
	double	dthresh,dthresh2, t1t2, darea,
	mx,my, mx2,my2,mxy,  tv,
	xm,ym, xm2,ym2,xym,
	temp,temp2, theta,pmx2,pmy2;
	int		x, y, xmin,xmax, ymin,ymax,area2, fdnpix, dnpix;

	float rv;
	int pixIndex;

	dthresh = d_cuttedDthreshArray[tid];

	if (plistexist_dthresh) //false
		minthresh = BIG;
	else	//true
		minthresh = 0.0;

	fdnpix = dnpix = 0;
	rv = 0.0;
	peak = cpeak = -BIG;
	ymin = xmin = 2*MAXPICSIZE;    /* to be really sure!! */
	ymax = xmax = 0;

	for(int i=index; i<index+pixcount; i++)
	{
		pixIndex = pixIndexArray[i];
		x = pixIndex % width;
		y = pixIndex / width;

		d_labelArray[pixIndex] = tid;

		if(peak < (val = d_pixelArray[pixIndex]))
			peak = val;
		if(cpeak < (cval = d_cdPixArray[pixIndex]))
			cpeak = cval;

		rv += cval;
		if (xmin > x)
			xmin = x;
		if (xmax < x)
			xmax = x;
		if (ymin > y)
			ymin = y;
		if (ymax < y)
			ymax = y;
		fdnpix++;
	}

	if(level > 0)
		flag |= OBJ_MERGED;
	if(xmin == 0 || xmax == width-1 || ymin == 0 || ymax == height-1)
		flag |= OBJ_TRUNC;

	/* copy some data to "obj" structure */
	//obj->fdnpix = (LONG)fdnpix;
	d_fdflux[tid] = rv;
	d_fdpeak[tid] = cpeak;
	d_dpeak[tid]  = peak;
	d_xmin[tid] = xmin;
	d_xmax[tid] = xmax;
	d_ymin[tid] = ymin;
	d_ymax[tid] = ymax;

	if (analyse_type & ANALYSE_FULL)
	{
		mx = my = tv = 0.0;
		mx2 = my2 = mxy = 0.0;
		dthresh2 = (dthresh + peak)/2.0;
		area2 = 0;

		for(int i=index; i<index+pixcount; i++)
		{
			pixIndex = pixIndexArray[i];
			x = pixIndex % width - xmin; /* avoid roundoff errors on big images */
			y = pixIndex / width - ymin; /* avoid roundoff errors on big images */

			cval = cval = d_cdPixArray[pixIndex];
			tv += (val = d_pixelArray[pixIndex]);

			if (val>dthresh)
				dnpix++;
			if (val > dthresh2)
				area2++;
			mx += cval * x;
			my += cval * y;
			mx2 += cval * x*x;
			my2 += cval * y*y;
			mxy += cval * x*y;
		}

		/*----- compute object's properties */
		xm = mx / rv;			/* mean x */
		ym = my / rv;			/* mean y */

		/*-- In case of blending, use previous barycenters (why?)*/
		xm2 = mx2 / rv - xm * xm;	/* variance of x */
		ym2 = my2 / rv - ym * ym;	/* variance of y */
		xym = mxy / rv - xm * ym;	/* covariance */

		/* Handle fully correlated x/y (which cause a singularity...) */
		if ((temp2=xm2*ym2-xym*xym)<0.00694)
		{
			xm2 += 0.0833333;
			ym2 += 0.0833333;
			temp2 = xm2*ym2-xym*xym;
			d_singuflag[tid] = 1;
		}
		else;
		d_singuflag[tid] = 0;

		if ((fabs(temp=xm2-ym2)) > 0.0)
			theta = atan2(2.0 * xym,temp) / 2.0;
		else
			theta = PI/4.0;

		temp = sqrt(0.25*temp*temp+xym*xym);
		pmy2 = pmx2 = 0.5*(xm2+ym2);
		pmx2+=temp;
		pmy2-=temp;

		d_dnpix[tid] = (flag & OBJ_OVERFLOW)? pixcount:dnpix;
		d_dflux[tid] = tv;
		d_mx[tid] = xm+xmin;	/* add back xmin */
		d_my[tid] = ym+ymin;	/* add back ymin */
		d_mx2[tid] = xm2;
		d_my2[tid] = ym2;
		d_mxy[tid] = xym;
		d_a[tid] = (float)sqrt(pmx2);
		d_b[tid] = (float)sqrt(pmy2);
		d_theta[tid] = theta*180.0/PI;

		d_cxx[tid] = (float)(ym2/temp2);
		d_cyy[tid] = (float)(xm2/temp2);
		d_cxy[tid] = (float)(-2*xym/temp2);

		//if(tid == 0)
		//	printf("cxx,cyy, cxy is %f, %f, %f", d_cxx[tid], d_cyy[tid], d_cxy[tid]);

		darea = (double)area2 - dnpix;
		t1t2 = dthresh/dthresh2;
		if (t1t2>0.0 /* && !prefs.dweight_flag */)
		{
			d_abcor[tid] = (darea<0.0?darea:-1.0)/(2*PI*log(t1t2<1.0?t1t2:0.99)*d_a[tid]*d_b[tid]);
			if (d_abcor[tid] > 1.0)
				d_abcor[tid] = 1.0;
		}
		else
			d_abcor[tid] = 1.0;
	}
	d_cuttedObjFlagArray[tid] = flag;

	float dist = pixcount/(2*PI*d_abcor[tid]*d_a[tid]*d_b[tid]);
	d_amp[tid] = dist<70.0? thresh*exp(dist) : 4.0*d_fdpeak[tid];

	/* ------------ limitate expansion ! */
	if (d_amp[tid]>4.0*d_fdpeak[tid])
		d_amp[tid] = 4.0*d_fdpeak[tid];
}

__global__ void preanalyse_robust_kernel(
		unsigned int *d_cuttedObjLevelArray,
		unsigned int *d_cuttedObjIndexArray,
		unsigned int *d_cuttedPixCountArray,
		unsigned int *d_cuttedObjFlagArray,
		float		 *d_cuttedDthreshArray,
		unsigned int *d_finalPixelIndexArray,
		float		 *d_cdPixArray,
		float		 *d_pixelArray,
		/* output starts */
		unsigned int *d_xmin,
		unsigned int *d_xmax,
		unsigned int *d_ymin,
		unsigned int *d_ymax,
		unsigned int *d_dnpix,
		double *d_mx,
		double *d_my,
		double *d_mx2,
		double *d_my2,
		double *d_mxy,
		float *d_cxx,
		float *d_cxy,
		float *d_cyy,
		float *d_a,
		float *d_b,
		float *d_theta,
		float *d_abcor,
		float *d_fdpeak,
		float *d_dpeak,
		float *d_fdflux,
		float *d_dflux,
		char  *d_singuflag,
		/* output ends */
		int width,
		int height,
		float thresh,
		int plistexist_dthresh,
		int analyse_type,
		unsigned int numobj)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= numobj)
		return;

	//get the level the object belongs to.
	unsigned int level = d_cuttedObjLevelArray[tid];

	if(level == 0)
		return;

	//get the index of the object in original pix index array.
	unsigned int index = d_cuttedObjIndexArray[tid];
	//get the number of pixels of the object
	unsigned int pixcount = d_cuttedPixCountArray[tid];

	unsigned int flag = d_cuttedObjFlagArray[tid];

	//start to compute the attributes
	/*-----  initialize stacks and bounds */
	float	peak, cpeak, val, cval, minthresh, thresht;
	double	dthresh,dthresh2, t1t2, darea,
	mx,my, mx2,my2,mxy,  tv,
	xm,ym, xm2,ym2,xym,
	temp,temp2, theta,pmx2,pmy2;
	int		x, y, xmin,xmax, ymin,ymax,area2, fdnpix, dnpix;

	float rv;
	int pixIndex;

	dthresh = d_cuttedDthreshArray[tid];

	if (plistexist_dthresh) //false
		minthresh = BIG;
	else	//true
		minthresh = 0.0;

	fdnpix = dnpix = 0;
	rv = 0.0;
	peak = cpeak = -BIG;
	ymin = xmin = 2*MAXPICSIZE;    /* to be really sure!! */
	ymax = xmax = 0;

	for(int i=index; i<index+pixcount; i++)
	{
		pixIndex = d_finalPixelIndexArray[i];
		x = pixIndex % width;
		y = pixIndex / width;

		if(peak < (val = d_pixelArray[pixIndex]))
			peak = val;
		if(cpeak < (cval = d_cdPixArray[pixIndex]))
			cpeak = cval;

		rv += cval;
		if (xmin > x)
			xmin = x;
		if (xmax < x)
			xmax = x;
		if (ymin > y)
			ymin = y;
		if (ymax < y)
			ymax = y;
		fdnpix++;
	}

	//if(level > 0)
	//	flag |= OBJ_MERGED;
	if(xmin == 0 || xmax == width-1 || ymin == 0 || ymax == height-1)
		flag |= OBJ_TRUNC;

	/* copy some data to "obj" structure */
	//obj->fdnpix = (LONG)fdnpix;
	d_fdflux[tid] = rv;
	d_fdpeak[tid] = cpeak;
	d_dpeak[tid]  = peak;
	d_xmin[tid] = xmin;
	d_xmax[tid] = xmax;
	d_ymin[tid] = ymin;
	d_ymax[tid] = ymax;

	if (analyse_type & ANALYSE_FULL)
	{
		mx = my = tv = 0.0;
		mx2 = my2 = mxy = 0.0;
		dthresh2 = (dthresh + peak)/2.0;
		area2 = 0;

		for(int i=index; i<index+pixcount; i++)
		{
			pixIndex = d_finalPixelIndexArray[i];
			x = pixIndex % width - xmin; /* avoid roundoff errors on big images */
			y = pixIndex / width - ymin; /* avoid roundoff errors on big images */

			cval = cval = d_cdPixArray[pixIndex];
			tv += (val = d_pixelArray[pixIndex]);

			if (val>dthresh)
				dnpix++;
			if (val > dthresh2)
				area2++;
			mx += cval * x;
			my += cval * y;
			mx2 += cval * x*x;
			my2 += cval * y*y;
			mxy += cval * x*y;
		}

		/*----- compute object's properties */
		xm = mx / rv;			/* mean x */
		ym = my / rv;			/* mean y */

		/*-- In case of blending, use previous barycenters (why?)*/
		if (flag & OBJ_MERGED)
		{
			double	xn,yn;

			xn = d_mx[tid]-xmin;
			yn = d_my[tid]-ymin;
			xm2 = mx2 / rv + xn*xn - 2*xm*xn;
			ym2 = my2 / rv + yn*yn - 2*ym*yn;
			xym = mxy / rv + xn*yn - xm*yn - xn*ym;
			xm = xn;
			ym = yn;
		}
		else
		{
			xm2 = mx2 / rv - xm * xm;	/* variance of x */
			ym2 = my2 / rv - ym * ym;	/* variance of y */
			xym = mxy / rv - xm * ym;	/* covariance */
		}

		/* Handle fully correlated x/y (which cause a singularity...) */
		if ((temp2=xm2*ym2-xym*xym)<0.00694)
		{
			xm2 += 0.0833333;
			ym2 += 0.0833333;
			temp2 = xm2*ym2-xym*xym;
			d_singuflag[tid] = 1;
		}
		else;
		d_singuflag[tid] = 0;

		if ((fabs(temp=xm2-ym2)) > 0.0)
			theta = atan2(2.0 * xym,temp) / 2.0;
		else
			theta = PI/4.0;

		temp = sqrt(0.25*temp*temp+xym*xym);
		pmy2 = pmx2 = 0.5*(xm2+ym2);
		pmx2+=temp;
		pmy2-=temp;

		d_dnpix[tid] = (flag & OBJ_OVERFLOW)? pixcount:dnpix;
		d_dflux[tid] = tv;
		d_mx[tid] = xm+xmin;	/* add back xmin */
		d_my[tid] = ym+ymin;	/* add back ymin */
		d_mx2[tid] = xm2;
		d_my2[tid] = ym2;
		d_mxy[tid] = xym;
		d_a[tid] = (float)sqrt(pmx2);
		d_b[tid] = (float)sqrt(pmy2);
		d_theta[tid] = theta*180.0/PI;

		d_cxx[tid] = (float)(ym2/temp2);
		d_cyy[tid] = (float)(xm2/temp2);
		d_cxy[tid] = (float)(-2*xym/temp2);

		darea = (double)area2 - dnpix;
		t1t2 = dthresh/dthresh2;
		if (t1t2>0.0 /* && !prefs.dweight_flag */)
		{
			d_abcor[tid] = (darea<0.0?darea:-1.0)/(2*PI*log(t1t2<1.0?t1t2:0.99)*d_a[tid]*d_b[tid]);
			if (d_abcor[tid] > 1.0)
				d_abcor[tid] = 1.0;
		}
		else
			d_abcor[tid] = 1.0;
	}
	d_cuttedObjFlagArray[tid] = flag;
}

__device__ float _back(float *d_back, int x, int y, int nx, int ny, int backw, int backh)
{
	int		xl,yl, pos;
	double	dx,dy, cdx;
	float	*b, b0,b1,b2,b3;

	b = d_back;

	dx = (double)x/backw - 0.5;
	dy = (double)y/backh - 0.5;
	dx -= (xl = (int)dx);
	dy -= (yl = (int)dy);

	if (xl<0)
	{
		xl = 0;
		dx -= 1.0;
	}
	else if (xl>=nx-1)
	{
		xl = nx<2 ? 0 : nx-2;
		dx += 1.0;
	}

	if (yl<0)
	{
		yl = 0;
		dy -= 1.0;
	}
	else if (yl>=ny-1)
	{
		yl = ny<2 ? 0 : ny-2;
		dy += 1.0;
	}

	pos = yl*nx + xl;
	cdx = 1 - dx;

	b0 = *(b+=pos);		/* consider when nbackx or nbacky = 1 */
	b1 = nx<2? b0:*(++b);
	b2 = ny<2? *b:*(b+=nx);
	b3 = nx<2? *b:*(--b);

	return (float)((1-dy)*(cdx*b0 + dx*b1) + dy*(dx*b2 + cdx*b3));
}

//prefs.pback_type == LOCAL case is not handled,
//we only handle prefs.pback_type == GLOBAL
__global__ void analyse_kernel(
		unsigned int 	*d_cuttedObjLabelArray,
		unsigned int 	*d_cuttedObjIndexArray,
		unsigned int	*d_cuttedPixCountArray,
		unsigned int	*d_finalPixelIndexArray,
		unsigned int	*d_npix,
		unsigned int 	*d_peakx,
		unsigned int 	*d_peaky,
		int			**d_iso,
		float		*d_back,
		float		*d_sigma,
		float		*d_cdPixArray,
		float		*d_pixelArray,
		float		*d_a,
		float		*d_b,
		float		*d_fdflux,
		float		*d_flux,
		float		*d_fluxerr,
		float		*d_peak,
		float 		*d_bkg,
		float 		*d_dbkg,
		float		*d_sigbkg,
		float		*d_mthresh,
		float		*d_thresh,
		float		*d_fwhm,
		char		*d_singuflag,
		double		*d_mx,
		double		*d_my,
		double		*d_poserr_mx2,
		double		*d_poserr_my2,
		double		*d_poserr_mxy,
		/* flag object vars */
		double		flagobj_poserr_mx2,
		int			flagobj_peakx,
		int			flagobj_iso_0,
		float		flagobj_fwhm,
		/* field global vars */
		float		backsig,
		double		gain,
		double		ngamma,
		float		thresh,
		float		dthresh,
		double		satur_level,
		/*	constant vars */
		int			plistexist_var,
		int			plistoff_var,
		int			prefs_weightgain_flag,
		int			prefs_clean_flag,
		int			detect_type,
		int 		nbackx,
		int 		nbacky,
		int			backw,
		int			backh,
		int			width,
		int 		numobj)
{
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= numobj)
		return;

	d_cuttedObjLabelArray[tid] = tid;

	d_bkg[tid] = _back(d_back,
			(int)(d_mx[tid]+0.5),
			(int)(d_my[tid]+0.5),
			nbackx, nbacky, backw, backh);

	d_dbkg[tid] = 0.0;

	//prefs.pback_type == GLOBAL
	d_sigbkg[tid] = backsig;

	///////////////////examineiso

	int			i,j,k,h, photoflag,area,errflag, cleanflag,
	pospeakflag, gainflag;
	double		tv,sigtv,
	esum, emx2,emy2,emxy, err,backnoise2,dbacknoise2,
	xm,ym, x,y,var,var2;
	float		*heapt,*heapj,*heapk, swap;
	PIXTYPE		pix, cdpix, tpix, peak,cdpeak;

	/*(?) is it ok to remove static */
	PIXTYPE		threshs[NISO];

	//we change dynamic allocate to static
	float 		heap[EXT_MINAREA];

	/* Prepare computation of positional error */
	esum = emx2 = emy2 = emxy = 0.0;
	if ((errflag=(*((char *)&flagobj_poserr_mx2))))
	{
		dbacknoise2 = backsig*backsig;
		xm = d_mx[tid];
		ym = d_my[tid];
	}
	else
		xm = ym = dbacknoise2 = 0.0;	/* to avoid gcc -Wall warnings */

	pospeakflag = (*((char *)&flagobj_peakx));
	//gain = field->gain;
	//ngamma = field->ngamma;
	photoflag = (detect_type==1/*PHOTO*/);

	gainflag = (plistexist_var) && prefs_weightgain_flag;

	h = EXT_MINAREA;

	/* Prepare selection of the heap selection for CLEANing */
	if ((cleanflag = prefs_clean_flag))
	{
		if (d_cuttedPixCountArray[tid] < EXT_MINAREA)
		{
			d_mthresh[tid] = 0.0;
			cleanflag = 0;
			//heapt = heap = NULL;		/* to avoid gcc -Wall warnings */
		}
		else
		{
			//heap = (float*)malloc(minarea*sizeof(float));
			heapt = heap;
		}
	}
	else
	{
		d_mthresh[tid] = 0.0;
		//heapt = heap = NULL;		/* to avoid gcc -Wall warnings */
	}

	/* Measure essential isophotal parameters in the measurement image... */
	tv = sigtv = 0.0;
	var = backnoise2 = backsig*backsig;
	peak = -BIG;
	cdpeak = -BIG;
	//thresh = field->thresh;
	//dthresh = dfield->dthresh;
	area = 0;

	unsigned int index = d_cuttedObjIndexArray[tid];
	unsigned int pixcount = d_cuttedPixCountArray[tid];
	int pixIndex;

	//added on gpu
	if(!pospeakflag)
	{
		d_peakx[tid] = 0;
		d_peaky[tid] = 0;
	}

	for(i=index; i<index+pixcount; i++)
	{
		pixIndex = d_finalPixelIndexArray[i];
		x = pixIndex % width;
		y = pixIndex / width;

		pix = d_pixelArray[pixIndex]; //PLIST(pixt,value);
		if (pix>peak)
			peak = pix;

		cdpix=d_cdPixArray[pixIndex]; //PLISTPIX(pixt,cdvalue);

		if (pospeakflag && cdpix>cdpeak)
		{
			cdpeak=cdpix;
			d_peakx[tid] =  x + 1;
			d_peaky[tid] =  y + 1;
		}
		//if (plistexist_var)
		//	var = (*((PIXTYPE *)((pixt)+plistoff_var)));//PLISTPIX(pixt, var);
		if (photoflag)
		{
			pix = exp(pix/ngamma);
			var2 = pix*pix*var;
		}
		else
			var2 = var;

		if (gainflag && pix>0.0 && gain>0.0)
			var2 += pix/gain*var/backnoise2;

		sigtv += var2;

		if (pix>thresh)
			area++;
		tv += pix;

		if (errflag)
		{
			err = dbacknoise2;
			if (gain>0.0 && cdpix>0.0)
				err += cdpix/gain;
			x = x - xm;
			y = y - ym;
			esum += err;
			emx2 += err*x*x;
			emy2 += err*y*y;
			emxy += err*x*y;
		}

		/*-- Find the minareath pixel in decreasing intensity for CLEANing */
		if (cleanflag)
		{
			tpix = cdpix - (/*PLISTEXIST(dthresh)?PLISTPIX(pixt, dthresh):*/dthresh);

//			if (h>0)
//				*(heapt++) = (float)tpix;
//			else if (h)
//			{
//				if ((float)tpix>*heap)
//				{
//					*heap = (float)tpix;
//					for (j=0; (k=(j+1)<<1)<=EXT_MINAREA; j=k)
//					{
//						heapk = heap+k;
//						heapj = heap+j;
//						if (k != EXT_MINAREA && *(heapk-1) > *heapk)
//						{
//							heapk++;
//							k++;
//						}
//						if (*heapj <= *(--heapk))
//							break;
//						swap = *heapk;
//						*heapk = *heapj;
//						*heapj = swap;
//					}
//				}
//			}
//			else
//				hmedian1(heap, EXT_MINAREA);
//			h--;

			if (h>0)
				*(heapt++) = (float)tpix;
			else
			{
				hmedian1(heap, EXT_MINAREA);
				if((float)tpix>*heap)
					*heap = (float)tpix;
			}
			h--;
		}
	}

	/* Flagging from the flag-map */
	//if (PLISTEXIST(flag)) //false
	//	getflags(obj, pixel);

	if (cleanflag)
	{
		hmedian1(heap, EXT_MINAREA);
		d_mthresh[tid] = *heap;
		//free(heap);
	}

	if (errflag)
	{
		double	flux2;

		flux2 = d_fdflux[tid]*d_fdflux[tid];
		/*-- Estimation of the error ellipse moments: we re-use previous variables */
		emx2 /= flux2;	/* variance of xm */
		emy2 /= flux2;	/* variance of ym */
		emxy /= flux2;	/* covariance */

		/*-- Handle fully correlated profiles (which cause a singularity...) */
		esum *= 0.08333/flux2;
		if (d_singuflag[tid] && (emx2*emy2-emxy*emxy) < esum*esum)
		{
			emx2 += esum;
			emy2 += esum;
		}

		d_poserr_mx2[tid] = emx2;
		d_poserr_my2[tid] = emy2;
		d_poserr_mxy[tid] = emxy;
	}

	if (photoflag)
	{
		tv = ngamma*(tv-d_cuttedPixCountArray[tid]*exp(d_dbkg[tid]/ngamma));
		sigtv /= ngamma*ngamma;
	}
	else
	{
		tv -= d_cuttedPixCountArray[tid]*d_dbkg[tid];
		if (!gainflag && gain > 0.0 && tv>0.0)
			sigtv += tv/gain;
	}

	d_npix[tid] = area;
	d_flux[tid] = tv;
	d_fluxerr[tid] = sigtv;
	//d_peak[tid] = peak;

	d_thresh[tid] = thresh - d_dbkg[tid];
	d_peak[tid] = peak - d_dbkg[tid];

	/* Initialize isophotal thresholds so as to sample optimally the full profile*/

	if (*((char *)&flagobj_iso_0))
	{
		int	**iso;
		PIXTYPE	*thresht;

		//memset(obj->iso, 0, NISO*sizeof(int));
		if (detect_type == 1/*PHOTO*/)
			for (i=0; i<NISO; i++)
				threshs[i] = d_thresh[tid] + (d_peak[tid]-d_thresh[tid])*i/NISO;
		else
		{
			if (d_peak[tid]>0.0 && d_thresh[tid]>0.0)
				for (i=0; i<NISO; i++)
					threshs[i] = d_thresh[tid]*pow((double)d_peak[tid]/d_thresh[tid], (double)i/NISO);
			else
				for (i=0; i<NISO; i++)
					threshs[i] = 0.0;
		}

		for(j=index; j<index+pixcount; j++)
		{
			pixIndex = d_finalPixelIndexArray[j];
			//pix = d_pixelArray[pixIndex]; //PLIST(pixt,value);

			for (i=NISO,iso=d_iso,thresht=threshs;
					i-- && d_pixelArray[pixIndex]>*(thresht++);)
				((*(iso++))[tid])++;
		}
	}

	/* Put objects in "segmentation check-image" */
	/*
	if ((check = prefs.check[CHECK_SEGMENTATION]))
		for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
			((ULONG *)check->pix)[check->width*PLIST(pixt,y)+PLIST(pixt,x)]
			                      = (ULONG)obj->number;

	if ((check = prefs.check[CHECK_OBJECTS]))
		for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
			((PIXTYPE *)check->pix)[check->width*PLIST(pixt,y)+PLIST(pixt,x)]
			                        = PLIST(pixt,value);
	 */


	/* Compute the FWHM of the object */
	if (*((char *)&flagobj_fwhm))
	{
		PIXTYPE	thresh0;

		thresh0 = d_peak[tid]/5.0;
		if (thresh0<d_thresh[tid])
			thresh0 = d_thresh[tid];
		if (thresh0>0.0)
		{
			double	mx,my, s,sx,sy,sxx,sxy, dx,dy,d2, lpix,pix, b, inverr2, sat,
			dbkg, d, bmax;

			mx = d_mx[tid];
			my = d_my[tid];
			dbkg = d_dbkg[tid];
			sat = (double)(satur_level - dbkg);
			s = sx = sy = sxx = sxy = 0.0;

			for(j=index; j<index+pixcount; j++)
			{
				pixIndex = d_finalPixelIndexArray[j];
				x = pixIndex % width;
				y = pixIndex / width;

				pix = d_pixelArray[pixIndex] - dbkg; //PLIST(pixt,value);

				if (pix>thresh0 && pix<sat)
				{
					dx = x - mx;
					dy = y - my;
					lpix = log(pix);
					inverr2 = pix*pix;
					s += inverr2;
					d2 = dx*dx+dy*dy;
					sx += d2*inverr2;
					sxx += d2*d2*inverr2;
					sy += lpix*inverr2;
					sxy += lpix*d2*inverr2;
				}
			}
			d = s*sxx-sx*sx;

			if (fabs(d)>0.0)
			{
				b = -(s*sxy-sx*sy)/d;
				if (b<(bmax = 1/(13*d_a[tid]*d_b[tid])))	/* to have FWHM # 6 sigma */
					b = bmax;
				d_fwhm[tid] = (float)(1.6651/sqrt(b));
				/*----- correction for undersampling effects (established from simulations) */
				if (d_fwhm[tid] > 0.0)
					d_fwhm[tid] -= 1/(4*d_fwhm[tid]);
			}
			else
				d_fwhm[tid] = 0.0;

		}
		else
			d_fwhm[tid] = 0.0;
	}
}


__global__ void removepixel_kernel(float *d_pixelArray,
		float *d_cdPixArray,
		float threshold,
		int width,
		int height)
{
	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	if(tid_x < width && tid_y < height) {

		const int tid = tid_y * width + tid_x;

		if(d_cdPixArray[tid] >= threshold)
		{
			d_pixelArray[tid] = -BIG;
		}
	}
}

__global__ void endobject_kernel(

		float		 *d_pixelArray,
		int			 *d_labelArray, //allocation map
		unsigned int *d_fdnpix,
		unsigned int *d_xmin,
		unsigned int *d_xmax,
		unsigned int *d_ymin,
		unsigned int *d_ymax,
		unsigned int *d_flag,

		int			**d_iso,
		float		*d_dthresh,
		float		*d_a,
		float		*d_b,
		float		*d_cxx,
		float		*d_cyy,
		float		*d_cxy,
		float		*d_peak,
		float		*d_bkg,
		float 		*d_dbkg,
		float		*d_fwhm,
		float		*d_theta,
		double		*d_mx,
		double		*d_my,
		double		*d_mxw,
		double		*d_myw,
		double		*d_alphas,
		double		*d_deltas,
		double		*d_alpha2000,
		double		*d_delta2000,
		double		*d_poserr_mx2,
		double		*d_poserr_my2,
		double		*d_poserr_mxy,
		double 		*d_posx,
		double		*d_posy,
		float		*d_elong,
		float		*d_ellip,
		float		*d_polar,
		float		*d_sprob,
		float		*d_sposx,
		float		*d_sposy,
		//float		*d_poserr_a,
		//float		*d_poserr_b,
		//float		*d_poserr_theta,
		//float		*d_poserr_cxx,
		//float		*d_poserr_cyy,
		//float		*d_poserr_cxy,
		float		*d_flux,
		float		*d_flux_iso,
		float		*d_fluxerr,
		float		*d_fluxerr_iso,
		float		*d_flux_isocor,
		float		*d_fluxerr_isocor,
		float		*d_mag_iso,
		float		*d_magerr_iso,
		float		*d_kronfactor,
		float		*d_flux_auto,
		float		*d_fluxerr_auto,
		float		*d_mag_auto,
		float		*d_magerr_auto,
		float		**d_flux_aper,
		float		**d_fluxerr_aper,
		float		**d_mag_aper,
		float		**d_magerr_aper,

		float 		flagobj2_mag_iso,
		float 		flagobj2_magerr_iso,
		float		flagobj2_poserr_a,
		float		flagobj2_poserr_cxx,
		float		flagobj2_elong,
		float		flagobj2_ellip,
		float		flagobj2_polar,
		float		flagobj2_flux_isocor,
		float		flagobj2_fluxerr_isocor,
		float		flagobj2_flux_auto,
		float		flagobj2_flux_petro,
		float       flagobj2_sprob,
		float		flagobj2_mag_auto,
		float		flagobj2_magerr_auto,
		float		*flagobj2_flux_aper,
		double		flagobj2_alpha2000,
		double		flagobj2_mxw,
		prefstruct		*d_prefs,
		brainstruct 	*d_brain,
		outobjstruct	*d_objlist,
		wcsstruct		*d_wcs,
		/*---field vars---*/
		int			fymin,
		int			fymax,
		int			width,
		int			height,
		float		backsig,
		float		thresh,
		double		satur_level,
		double		pixscale,
		double		ngamma,
		double		gain,
		int numobj)
{

	const int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if(tid >= numobj)
		return;

	////////////////

	int		i,j, ix,iy,selecflag, newnumber,nsub;

	/* Source position */
	d_sposx[tid] = (float)(d_posx[tid] = d_mx[tid]+1.0); /* That's standard FITS */
	d_sposy[tid] = (float)(d_posy[tid] = d_my[tid]+1.0);

	/* Integer coordinates */
	ix=(int)(d_mx[tid]+0.49999);
	iy=(int)(d_my[tid]+0.49999);

	/* Association */
	//if (prefs.assoc_flag) //false
	//	obj2->assoc_number = do_assoc(field, obj2->sposx, obj2->sposy);

	//if (prefs.assoc_flag && prefs.assocselec_type!=ASSOCSELEC_ALL) //false
	//	selecflag = (prefs.assocselec_type==ASSOCSELEC_MATCHED)?
	//			obj2->assoc_number:!obj2->assoc_number;
	//else
	selecflag = 1;

	if (selecflag)
	{
//		/*-- Paste back to the image the object's pixels if BLANKing is on */
//		if (d_prefs->blank_flag)
//		{
//			pasteimage(field, obj->blank, obj->subw, obj->subh,
//					obj->subx, obj->suby);
//			if (obj->dblank)
//				pasteimage(dfield, obj->dblank, obj->subw, obj->subh,
//						obj->subx, obj->suby);
//		}

		/*------------------------- Error ellipse parameters ------------------------*/
//		if (*((char *)&flagobj2_poserr_a))	//false
//		{
//			double	pmx2,pmy2,temp,theta;
//
//			if (fabs(temp=d_poserr_mx2[tid]-d_poserr_my2[tid]) > 0.0)
//				theta = atan2(2.0 * d_poserr_mxy[tid],temp) / 2.0;
//			else
//				theta = PI/4.0;
//
//			temp = sqrt(0.25*temp*temp+d_poserr_mxy[tid]*d_poserr_mxy[tid]);
//			pmy2 = pmx2 = 0.5*(d_poserr_mx2[tid]+d_poserr_my2[tid]);
//			pmx2+=temp;
//			pmy2-=temp;
//
//			d_poserr_a[tid] = (float)sqrt(pmx2);
//			d_poserr_b[tid] = (float)sqrt(pmy2);
//			d_poserr_theta[tid] = theta*180.0/PI;
//		}

//		if (*((char *)&flagobj2_poserr_cxx))	//false
//		{
//			double	xm2,ym2, xym, temp;
//
//			xm2 = d_poserr_mx2[tid];
//			ym2 = d_poserr_my2[tid];
//			xym = d_poserr_mxy[tid];
//			d_poserr_cxx[tid] = (float)(ym2/(temp=xm2*ym2-xym*xym));
//			d_poserr_cyy[tid] = (float)(xm2/temp);
//			d_poserr_cxy[tid] = (float)(-2*xym/temp);
//		}

		/* ---- Aspect ratio */

		if (*((char *)&flagobj2_elong))
			d_elong[tid] = d_a[tid]/d_b[tid];

		if (*((char *)&flagobj2_ellip))
			d_ellip[tid] = 1-d_b[tid]/d_a[tid];

		if (*((char *)&flagobj2_polar))
			d_polar[tid] = (d_a[tid]*d_a[tid] - d_b[tid]*d_b[tid])
			/ (d_a[tid]*d_a[tid] + d_b[tid]*d_b[tid]);

		/*------------------------------- Photometry -------------------------------*/

		/*-- Convert the father of photom. error estimates from variance to RMS */
		d_flux_iso[tid] = d_flux[tid];
		d_fluxerr_iso[tid] = sqrt(d_fluxerr[tid]);

		if (*((char *)&flagobj2_flux_isocor))
		{
			computeisocorflux(
					d_fdnpix,
					d_dthresh,
					d_flux,
					d_fluxerr,
					d_flux_isocor,
					d_fluxerr_isocor,
					flagobj2_fluxerr_isocor,
					tid);
		}

		if (*((char *)&flagobj2_flux_aper))
			for (i=0; i<NAPER; i++)
				computeaperflux(
				d_pixelArray,
				d_labelArray, //allocation map
				d_flag,
				d_dbkg,
				d_flux_aper,
				d_fluxerr_aper,
				d_mag_aper,
				d_magerr_aper,
				d_mx,
				d_my,
				d_prefs,
				width,
				height,
				fymin,
				fymax,
				ngamma,
				gain,
				backsig,
				tid,
				i);

		if (*((char *)&flagobj2_flux_auto))
			computeautoflux(
					d_pixelArray,
					d_labelArray, //allocation map
					d_flag,
					d_a,
					d_b,
					d_cxx,
					d_cyy,
					d_cxy,
					d_dbkg,
					d_kronfactor,
					d_flux_auto,
					d_fluxerr_auto,
					d_mag_auto,
					d_magerr_auto,
					d_mx,
					d_my,
					d_prefs,
					width,
					height,
					fymin,
					fymax,
					ngamma,
					gain,
					backsig,
					flagobj2_mag_auto,
					flagobj2_magerr_auto,
					tid);


		if (*((char *)&flagobj2_flux_petro))	//false
			//computepetroflux(field, dfield, wfield, dwfield, obj);

		/*-- Growth curve */
		if (d_prefs->growth_flag)	//false
			//makeavergrowth(field, wfield, obj);

		/*--------------------------- Windowed barycenter --------------------------*/
		//if (*((char *)&flagobj2_winpos_x))	//false
			//compute_winpos(field, wfield, obj);

		/*-- What about the peak of the profile? */
		if (d_peak[tid]+d_bkg[tid] >= satur_level)
			d_flag[tid] |= OBJ_SATUR;

		/*-- Check-image CHECK_APERTURES option */

//		if ((check = prefs.check[CHECK_APERTURES])) //false
//		{
//			if (FLAG(obj2.flux_aper))
//				for (i=0; i<prefs.naper; i++)
//					sexcircle(check->pix, check->width, check->height,
//							obj->mx, obj->my, prefs.apert[i]/2.0, check->overlay);
//
//			if (FLAG(obj2.flux_auto))
//				sexellips(check->pix, check->width, check->height,
//						obj->mx, obj->my, obj->a*obj2->kronfactor,
//						obj->b*obj2->kronfactor, obj->theta,
//						check->overlay, obj->flag&OBJ_CROWDED);
//
//			if (FLAG(obj2.flux_petro))
//				sexellips(check->pix, check->width, check->height,
//						obj->mx, obj->my, obj->a*obj2->petrofactor,
//						obj->b*obj2->petrofactor, obj->theta,
//						check->overlay, obj->flag&OBJ_CROWDED);
//		}

		/* ---- Star/Galaxy classification */

		if (*((char *)&flagobj2_sprob))
		{
			double	fac2, input[10], output, fwhm;

			fwhm = d_prefs->seeing_fwhm;

			fac2 = fwhm/pixscale;
			fac2 *= fac2;

			input[j=0] = log10(d_iso[0][tid]? d_iso[0][tid]/fac2: 0.01);

			input[++j] = thresh>0.0?
					log10(d_peak[tid]>0.0? d_peak[tid]/thresh: 0.1) :-1.0;
			for (i=1; i<NISO; i++)
				input[++j] = log10(d_iso[i][tid]? d_iso[i][tid]/fac2: 0.01);
			input[++j] = log10(fwhm);

			neurresp(d_brain, input, &output);
			d_sprob[tid] = (float)output;
		}

		/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
		/*-- Put here your calls to "BLIND" custom functions. Ex:

    compute_myotherparams(obj);

--*/

		/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

//		newnumber = ++thecat.ntotal;
//		/*-- update segmentation map */
//		if ((check=prefs.check[CHECK_SEGMENTATION])) //false
//		{
//			ULONG	*pix;
//			ULONG	newsnumber = newnumber,
//					oldsnumber = obj->number;
//			int	dx,dx0,dy,dpix;
//
//			pix = (ULONG *)check->pix + check->width*obj->ymin + obj->xmin;
//			dx0 = obj->xmax-obj->xmin+1;
//			dpix = check->width-dx0;
//			for (dy=obj->ymax-obj->ymin+1; dy--; pix += dpix)
//				for (dx=dx0; dx--; pix++)
//					if (*pix==oldsnumber)
//						*pix = newsnumber;
//		}
//		obj->number = newnumber;

//		/*-- SOM fitting */
//		if (prefs.somfit_flag) //false
//		{
//			float    *input;
//
//			input = thesom->input;
//			copyimage(field,input,thesom->inputsize[0],thesom->inputsize[1],ix,iy);
//
//			if (thesom->nextrainput)
//			{
//				input += thesom->ninput-thesom->nextrainput;
//				*(input) = (obj->mx+1)/field->width;
//				*(input+1) = (obj->my+1)/field->height;
//			}
//
//			som_phot(thesom, obj->bkg, field->backsig,
//					(float)field->gain, obj->mx-ix, obj->my-iy,
//					FLAG(obj2.vector_somfit)?outobj2.vector_somfit:NULL, -1.0);
//			obj2->stderr_somfit = thesom->stderror;
//			obj2->flux_somfit = thesom->amp;
//			outobj2.fluxerr_somfit = thesom->sigamp;
//		}

//		if (FLAG(obj2.vignet)) //false
//			copyimage(field,outobj2.vignet,prefs.vignetsize[0],prefs.vignetsize[1],
//					ix,iy);
//
//		if (FLAG(obj2.vigshift)) //false
//			copyimage_center(field, outobj2.vigshift, prefs.vigshiftsize[0],
//					prefs.vigshiftsize[1], obj->mx, obj->my);

//		/*------------------------------- PSF fitting ------------------------------*/
//		nsub = 1;
//		if (prefs.psf_flag)	//false
//		{
//			if (prefs.dpsf_flag)
//				double_psf_fit(ppsf, field, wfield, obj, thepsf, dfield, dwfield);
//			else
//				psf_fit(thepsf, field, wfield, obj);
//			obj2->npsf = thepsfit->npsf;
//			if (prefs.psfdisplay_type == PSFDISPLAY_SPLIT)
//			{
//				nsub = thepsfit->npsf;
//				if (nsub<1)
//					nsub = 1;
//			}
//			else
//				for (j=0; j<thepsfit->npsf; j++)
//				{
//					if (FLAG(obj2.x_psf) && j<prefs.psf_xsize)
//						obj2->x_psf[j] = thepsfit->x[j];
//					if (FLAG(obj2.y_psf) && j<prefs.psf_ysize)
//						obj2->y_psf[j] = thepsfit->y[j];
//					if (FLAG(obj2.flux_psf) && j<prefs.psf_fluxsize)
//						obj2->flux_psf[j] = thepsfit->flux[j];
//					if (FLAG(obj2.magerr_psf) && j<prefs.psf_magerrsize)
//						obj2->magerr_psf[j] = obj2->fluxerr_psf[j]>0.0?
//								1.086*obj2->fluxerr_psf[j]/thepsfit->flux[j] : 99.0;
//					if (FLAG(obj2.fluxerr_psf) && j<prefs.psf_fluxerrsize)
//						obj2->fluxerr_psf[j] = obj2->fluxerr_psf[j];
//					if (FLAG(obj2.mag_psf) && j<prefs.psf_magsize)
//						obj2->mag_psf[j] = thepsfit->flux[j]>0.0?
//								prefs.mag_zeropoint -2.5*log10(thepsfit->flux[j]) : 99.0;
//				}
//		}

		/*----------------------------- Profile fitting -----------------------------*/
//		nsub = 1;
//		if (prefs.prof_flag) //false
//			profit_fit(theprofit, field, wfield, obj, obj2);

		/*--- Express everything in magnitude units */
		//computemags(field, obj);
		computemags(d_mag_iso,
				d_magerr_iso,
				d_flux_iso,
				d_fluxerr_iso,
				flagobj2_mag_iso,
				flagobj2_magerr_iso,
				d_prefs->mag_zeropoint,
				tid);

		/*-------------------------------- Astrometry ------------------------------*/
		if (d_prefs->world_flag) //true
			;
			/*computeastrom1(
				d_posx,
				d_posy,
				d_mxw,
				d_myw,
				d_alphas,
				d_deltas,
				d_alpha2000,
				d_delta2000,
				d_wcs,
				flagobj2_alpha2000,
				flagobj2_mxw,
				tid);*/

		/*-- Edit min and max coordinates to follow the FITS conventions */
		d_xmin[tid] += 1;
		d_ymin[tid] += 1;
		d_xmax[tid] += 1;
		d_ymax[tid] += 1;

		/*-- Go through each newly identified component */
		for (j=0; j<nsub; j++)
		{
//			if (prefs_psf_flag && prefs_psfdisplay_type == PSFDISPLAY_SPLIT)
//			{
//
//				if (*((char *)&flagobj2_x_psf))
//					obj2->x_psf[0] = thepsfit->x[j];
//				if (FLAG(obj2.y_psf))
//					obj2->y_psf[0] = thepsfit->y[j];
//				if (FLAG(obj2.flux_psf))
//					obj2->flux_psf[0] = thepsfit->flux[j]>0.0? thepsfit->flux[j]:0.0; /*?*/
//							if (FLAG(obj2.mag_psf))
//								obj2->mag_psf[0] = thepsfit->flux[j]>0.0?
//										prefs.mag_zeropoint -2.5*log10(thepsfit->flux[j]) : 99.0;
//							if (FLAG(obj2.magerr_psf))
//								obj2->magerr_psf[0]=
//										(thepsfit->flux[j]>0.0 && obj2->fluxerr_psf[j]>0.0) ? /*?*/
//												1.086*obj2->fluxerr_psf[j]/thepsfit->flux[j] : 99.0;
//							if (FLAG(obj2.fluxerr_psf))
//								obj2->fluxerr_psf[0]= obj2->fluxerr_psf[j];
//							if (j)
//								obj->number = ++thecat.ntotal;
//			}


		}
	}
	else
	{
//		/*-- Treatment of discarded detections */
//		/*-- update segmentation map */
//		if ((check=prefs.check[CHECK_SEGMENTATION]))
//		{
//			ULONG	*pix;
//			ULONG	oldsnumber = obj->number;
//			int	dx,dx0,dy,dpix;
//
//			pix = (ULONG *)check->pix + check->width*obj->ymin + obj->xmin;
//			dx0 = obj->xmax-obj->xmin+1;
//			dpix = check->width-dx0;
//			for (dy=obj->ymax-obj->ymin+1; dy--; pix += dpix)
//				for (dx=dx0; dx--; pix++)
//					if (*pix==oldsnumber)
//						*pix = 0;
//		}
	}

	/* Remove again from the image the object's pixels if BLANKing is on ... */
	/*-- ... and free memory */

//	if (prefs.blank_flag && obj->blank)
//	{
//		if (selecflag)
//		{
//			if (prefs.somfit_flag && (check=prefs.check[CHECK_MAPSOM]))
//				blankcheck(check, obj->blank, obj->subw, obj->subh,
//						obj->subx, obj->suby, (PIXTYPE)*(obj2->vector_somfit));
//
//		}
//		blankimage(field, obj->blank, obj->subw, obj->subh,
//				obj->subx, obj->suby, -BIG);
//		free(obj->blank);
//		if (obj->dblank)
//		{
//			blankimage(dfield, obj->dblank, obj->subw, obj->subh,
//					obj->subx, obj->suby, -BIG);
//			free(obj->dblank);
//		}
//	}

	(d_objlist + tid)->number = tid + 1;
	(d_objlist + tid)->alpha2000 = 0.0;		//ALPHA_J2000
	(d_objlist + tid)->delta2000 = 0.0;		//DELTA_J2000

	(d_objlist + tid)->sposx = d_sposx[tid];
	(d_objlist + tid)->sposy = d_sposy[tid];
	(d_objlist + tid)->flux_iso = d_flux_iso[tid];
	(d_objlist + tid)->fluxerr_iso = d_fluxerr_iso[tid];
	(d_objlist + tid)->mag_iso = d_mag_iso[tid];
	(d_objlist + tid)->magerr_iso = d_magerr_iso[tid];
	(d_objlist + tid)->flux_isocor = d_flux_isocor[tid];
	(d_objlist + tid)->fluxerr_isocor = d_fluxerr_isocor[tid];

	for (i=0; i<NAPER; i++)
	{
		(d_objlist + tid)->flux_aper[i] 	= d_flux_aper[i][tid];
		(d_objlist + tid)->fluxerr_aper[i] = d_fluxerr_aper[i][tid];
		(d_objlist + tid)->mag_aper[i] 	= d_mag_aper[i][tid];
		(d_objlist + tid)->magerr_aper[i] 	= d_magerr_aper[i][tid];
	}

	(d_objlist + tid)->flux_auto = d_flux_auto[tid];
	(d_objlist + tid)->fluxerr_auto = d_fluxerr_auto[tid];
	(d_objlist + tid)->mag_auto = d_mag_auto[tid];
	(d_objlist + tid)->magerr_auto = d_magerr_auto[tid];
	(d_objlist + tid)->theta = d_theta[tid];

	(d_objlist + tid)->poserr_mx2 = d_poserr_mx2[tid]; 	//ERRX2_IMAGE
	(d_objlist + tid)->poserr_my2 = d_poserr_my2[tid]; 	//ERRY2_IMAGE
	(d_objlist + tid)->poserr_mxy = d_poserr_mxy[tid]; 	//ERRXY_IMAGE
	(d_objlist + tid)->flag = d_flag[tid];				//FLAGS
	(d_objlist + tid)->elong = d_elong[tid];			//ELONGATION
	(d_objlist + tid)->ellip = d_ellip[tid];			//ELLIPTICITY
	(d_objlist + tid)->sprob = d_sprob[tid];			//CLASS_STAR
	(d_objlist + tid)->peak   = d_peak[tid];			//FLUX_MAX
	(d_objlist + tid)->kronfactor = d_kronfactor[tid]; 	//KRON_RADIUS
	(d_objlist + tid)->bkg = d_bkg[tid];				//BACKGROUND
	(d_objlist + tid)->fwhm = d_fwhm[tid];				//FWHM_IMAGE


	return;

}

//iso memset 0 before use
extern "C" void analyse_gpu(
		double		flagobj_poserr_mx2,
		int			flagobj_peakx,
		int			flagobj_iso_0,
		float		flagobj_fwhm,
		/* field global vars */
		float		backsig,
		double		gain,
		double		ngamma,
		float		thresh,
		float		dthresh,
		double		satur_level,
		/*	constant vars */
		int			plistexist_var,
		int			plistoff_var,
		int			prefs_weightgain_flag,
		int			prefs_ext_minarea,
		int			prefs_clean_flag,
		int			detect_type,
		int 		nbackx,
		int 		nbacky,
		int			backw,
		int			backh,
		int			width,
		int 		numobj)
{
	int grid = (numobj-1)/(MAX_THREADS_PER_BLK)+1;
	int block = (MAX_THREADS_PER_BLK);

	if(prefs_ext_minarea != EXT_MINAREA)
		printf("config and macro(EXT_MINAREA) do not match, please check!!!\n");

	int **p_iso;
	cudaMalloc((void**)&p_iso, NISO*sizeof(int*));
	cudaMemcpy(p_iso, d_iso, NISO*sizeof(int*), cudaMemcpyHostToDevice);

	float time;
	cudaEventRecord(start, 0);

	int plistexist_dthresh = 0;
	//run preanalyse again but add the new pixels in allocated map.
	preanalyse_robust_kernel<<<grid, block>>>(
			d_cuttedObjLevelArray,
			d_cuttedObjIndexArray,
			d_cuttedPixCountArray,
			d_cuttedObjFlagArray,
			d_cuttedDthreshArray,
			d_finalPixelIndexArray,
			d_cdPixArray,
			d_pixelArray,
			/* output starts */
			d_xmin,
			d_xmax,
			d_ymin,
			d_ymax,
			d_dnpix,
			d_mx,
			d_my,
			d_mx2,
			d_my2,
			d_mxy,
			d_cxx,
			d_cxy,
			d_cyy,
			d_a,
			d_b,
			d_theta,
			d_abcor,
			d_fdpeak,
			d_dpeak,
			d_fdflux,
			d_dflux,
			d_singuflag,
			/* output ends */
			width,
			height,
			thresh,
			plistexist_dthresh,
			ANALYSE_FULL|ANALYSE_ROBUST,
			numobj);

	//d_iso copy and memset
	analyse_kernel<<<grid, block>>>(
			d_cuttedObjLabelArray,
			d_cuttedObjIndexArray,
			d_cuttedPixCountArray,
			d_finalPixelIndexArray,
			d_npix,
			d_peakx,
			d_peaky,
			p_iso,
			d_mean,
			d_sigma,
			d_cdPixArray,
			d_pixelArray,
			d_a,
			d_b,
			d_fdflux,
			d_flux,
			d_fluxerr,
			d_peak,
			d_bkg,
			d_dbkg,
			d_sigbkg,
			d_mthresh,
			d_thresh,
			d_fwhm,
			d_singuflag,
			d_mx,
			d_my,
			d_poserr_mx2,
			d_poserr_my2,
			d_poserr_mxy,
			/* flag object vars */
			flagobj_poserr_mx2,
			flagobj_peakx,
			flagobj_iso_0,
			flagobj_fwhm,
			/* field global vars */
			backsig,
			gain,
			ngamma,
			thresh,
			dthresh,
			satur_level,
			/*	constant vars */
			plistexist_var,
			plistoff_var,
			prefs_weightgain_flag,
			prefs_clean_flag,
			detect_type,
			nbackx,
			nbacky,
			backw,
			backh,
			width,
			numobj);

	cudaFree(p_iso);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

#ifdef DEBUG_CUDA
	printf("Time counsumed by analyse stage 1 is: %f\n", time);
#endif
}

__global__ void updateAllocateLabel(
		int 	*d_labelArray,
		float	*d_cdPixArray,
		unsigned int *d_masterIndex,
		float	thresh,
		int		width,
		int		height)
{
	const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

	if(tid_x < width && tid_y < height) {

		const int tid = tid_y * width + tid_x;
		if(d_labelArray[tid] > 0 && d_masterIndex[d_labelArray[tid]] != -1)
			d_labelArray[tid] = INT_MAX;

		if(d_cdPixArray[tid] >= thresh && d_labelArray[tid] == -1)
		{
			d_labelArray[tid] = INT_MAX;
		}
	}
}

extern "C" void endobject_gpu(
		prefstruct	*h_prefs,
		brainstruct *h_brain,
		obj2struct	*h_flagobj2,
		outobjstruct *h_outobj,
		wcsstruct	*h_wcs,
		/*---field vars---*/
		int			fymin,
		int			fymax,
		int			width,
		int			height,
		float		backsig,
		float		thresh,
		double		satur_level,
		double		pixscale,
		double		ngamma,
		double		gain,
		int 		numobj)
{
	if(h_prefs->naper != NAPER)
	{
		printf("config and macro(NAPER) definition do not match, please check!!!\n");
		exit(0);
	}

	float time;
	cudaEventRecord(start, 0);

	//update allocate map for pixels that are above threshold but don't belong to any object
	dim3 grid4((width - 1) / SQUARE_BLK_WIDTH + 1,
			(height - 1) / SQUARE_BLK_HEIGHT + 1);
	dim3 block4(SQUARE_BLK_WIDTH, SQUARE_BLK_HEIGHT);

	updateAllocateLabel<<<grid4, block4>>>(
			d_labelArray,
			d_cdPixArray,
			d_masterIndex,
			global_dthresh,
			width,
			height);

	prefstruct	*d_prefs;
	cudaMalloc((void**)&d_prefs, sizeof(prefstruct));
	cudaMemcpy(d_prefs, h_prefs, sizeof(prefstruct), cudaMemcpyHostToDevice);

	brainstruct	*d_brain;
	cudaMalloc((void**)&d_brain, sizeof(brainstruct));
	cudaMemcpy(d_brain, h_brain, sizeof(brainstruct), cudaMemcpyHostToDevice);

	float	**p_flux_aper, **p_fluxerr_aper, **p_mag_aper, **p_magerr_aper;

	cudaMalloc((void**)&p_flux_aper, 	NAPER*sizeof(float*));
	cudaMalloc((void**)&p_fluxerr_aper, NAPER*sizeof(float*));
	cudaMalloc((void**)&p_mag_aper, 	NAPER*sizeof(float*));
	cudaMalloc((void**)&p_magerr_aper, 	NAPER*sizeof(float*));
	cudaMemcpy(p_flux_aper, 	d_flux_aper, 	NAPER*sizeof(float*), cudaMemcpyHostToDevice);
	cudaMemcpy(p_fluxerr_aper, 	d_fluxerr_aper, NAPER*sizeof(float*), cudaMemcpyHostToDevice);
	cudaMemcpy(p_mag_aper, 		d_mag_aper, 	NAPER*sizeof(float*), cudaMemcpyHostToDevice);
	cudaMemcpy(p_magerr_aper, 	d_magerr_aper, 	NAPER*sizeof(float*), cudaMemcpyHostToDevice);

	int **p_iso;
	cudaMalloc((void**)&p_iso, NISO*sizeof(int*));
	cudaMemcpy(p_iso, d_iso, NISO*sizeof(int*), cudaMemcpyHostToDevice);

	outobjstruct *d_objlist;
	cudaMalloc((void**)&d_objlist, numobj*sizeof(outobjstruct));

	wcsstruct	*d_wcs;
	cudaMalloc((void**)&d_wcs, 	sizeof(wcsstruct));
	cudaMemcpy(d_wcs, h_wcs, sizeof(wcsstruct), cudaMemcpyHostToDevice);

	//wcsrev
	//   /* Initialize if required. */
	//   if (wcs->flag != WCSSET) {
	//      if (wcsset(lin->naxis, ctype, wcs)) return 1;
	//   }

	//linerev
	//   if (lin->flag != LINSET) {
	//      if (linset(lin)) return 1;
	//   }

	//celrev
	//	 if (cel->flag != CELSET) {
	//	 	if(celset(pcode, cel, prj)) return 1;
	//	}

	int grid = (numobj-1)/(MAX_THREADS_PER_BLK)+1;
	int block = (MAX_THREADS_PER_BLK);

	endobject_kernel<<<grid, block>>>(
			d_pixelArray,
			d_labelArray, //allocation map
			d_fdnpix,
			d_xmin,
			d_xmax,
			d_ymin,
			d_ymax,
			d_flag,
			p_iso,
			d_dthresh,
			d_a,
			d_b,
			d_cxx,
			d_cyy,
			d_cxy,
			d_peak,
			d_bkg,
			d_dbkg,
			d_fwhm,
			d_theta,
			d_mx,
			d_my,
			d_mxw,
			d_myw,
			d_alphas,
			d_deltas,
			d_alpha2000,
			d_delta2000,
			d_poserr_mx2,
			d_poserr_my2,
			d_poserr_mxy,
			d_posx,
			d_posy,
			d_elong,
			d_ellip,
			d_polar,
			d_sprob,
			d_sposx,
			d_sposy,
			d_flux,
			d_flux_iso,
			d_fluxerr,
			d_fluxerr_iso,
			d_flux_isocor,
			d_fluxerr_isocor,
			d_mag_iso,
			d_magerr_iso,
			d_kronfactor,
			d_flux_auto,
			d_fluxerr_auto,
			d_mag_auto,
			d_magerr_auto,
			p_flux_aper,
			p_fluxerr_aper,
			p_mag_aper,
			p_magerr_aper,
			h_flagobj2->mag_iso,
			h_flagobj2->magerr_iso,
			h_flagobj2->poserr_a,
			h_flagobj2->poserr_cxx,
			h_flagobj2->elong,
			h_flagobj2->ellip,
			h_flagobj2->polar,
			h_flagobj2->flux_isocor,
			h_flagobj2->fluxerr_isocor,
			h_flagobj2->flux_auto,
			h_flagobj2->flux_petro,
			h_flagobj2->sprob,
			h_flagobj2->mag_auto,
			h_flagobj2->magerr_auto,
			h_flagobj2->flux_aper,
			h_flagobj2->alpha2000,
			h_flagobj2->mxw,
			d_prefs,
			d_brain,
			d_objlist,
			d_wcs,
			/*---field vars---*/
			fymin,
			fymax,
			width,
			height,
			backsig,
			thresh,
			satur_level,
			pixscale,
			ngamma,
			gain,
			numobj);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

#ifdef DEBUG_CUDA
	printf("Time counsumed by cuda endobject is: %f\n", time);
#endif

	///////////////////////////////////
	time = 0.0;
	cudaEventRecord(start, 0);

	cudaMemcpy(h_outobj,  d_objlist,  numobj*sizeof(outobjstruct), cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

#ifdef DEBUG_CUDA
	printf("time consumed by cuda data transfer out is %f\n", time);
#endif

	cudaFree(d_prefs);
	cudaFree(d_brain);
	cudaFree(p_flux_aper);
	cudaFree(p_fluxerr_aper);
	cudaFree(p_mag_aper);
	cudaFree(p_magerr_aper);
	cudaFree(p_iso);
	cudaFree(d_objlist);

}

extern "C" void copyobjects(
		unsigned int *h_index,
		unsigned int *h_flag,
		unsigned int *h_fdnpix,
		float 		 *h_dthresh,

		unsigned int *h_xmin,
		unsigned int *h_xmax,
		unsigned int *h_ymin ,
		unsigned int *h_ymax,
		unsigned int *h_dnpix,
		unsigned int *h_npix,
		unsigned int *h_peakx,
		unsigned int *h_peaky ,

		float 	*h_bkg,
		float 	*h_dbkg ,
		float	*h_sigbkg,
		float 	*h_a,
		float 	*h_b,
		float 	*h_cxx,
		float 	*h_cxy,
		float 	*h_cyy,
		float 	*h_theta,
		float 	*h_abcor,
		float	*h_peak,
		float 	*h_dpeak,
		float 	*h_fdpeak ,
		float	*h_flux ,
		float 	*h_dflux,
		float 	*h_fdflux,
		float	*h_fluxerr,
		float	*h_thresh,
		float	*h_mthresh,
		float	*h_fwhm,

		double 	*h_mx,
		double 	*h_my,
		double 	*h_mx2,
		double 	*h_my2,
		double 	*h_mxy,
		double	*h_poserr_mx2,
		double	*h_poserr_my2,
		double	*h_poserr_mxy,
		char  	*h_singuflag,
		int		**h_iso,
		int total)
{
	/*
	objstruct *d_objectlist;
	cudaMalloc((void**)&d_objectlist, total*sizeof(objstruct));

	int grid = (total-1)/(MAX_THREADS_PER_BLK)+1;
	int block = (MAX_THREADS_PER_BLK);

	copyobj_kernel<<<grid, block>>>(
			d_index,
			d_flag,
			d_fdnpix,
			d_dthresh,
			d_xmin,
			d_xmax,
			d_ymin,
			d_ymax,
			d_dnpix,
			d_npix,
			d_peakx,
			d_peaky,

			d_bkg,
			d_dbkg,
			d_sigbkg,
			d_a,
			d_b,
			d_cxx,
			d_cxy,
			d_cyy,
			d_theta,
			d_abcor,
			d_peak,
			d_dpeak,
			d_fdpeak,
			d_flux,
			d_dflux,
			d_fdflux,
			d_fluxerr,
			d_thresh,
			d_mthresh,
			d_fwhm,
			d_mx,
			d_my,
			d_mx2,
			d_my2,
			d_mxy,
			d_poserr_mx2,
			d_poserr_my2,
			d_poserr_mxy,
			d_singuflag,
			d_iso,
			d_objectlist,
			total);

	cudaMemcpy(h_objectlist, d_objectlist, total*sizeof(objstruct), cudaMemcpyDeviceToHost);
	cudaFree(d_objectlist);
	*/

	cudaMemcpy(h_index,  d_index, total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_flag,   d_flag,  total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fdnpix, d_fdnpix,total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_dthresh, d_dthresh,total*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(h_xmin,	d_xmin,	total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_xmax,  d_xmax, total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_ymin, 	d_ymin,	total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_ymax,  d_ymax, total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_dnpix, d_dnpix,total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_npix, 	d_npix,	total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_peakx, d_peakx,total*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_peaky, d_peaky,total*sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaMemcpy(h_bkg, 	d_bkg,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_dbkg, 	d_dbkg,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_sigbkg,d_sigbkg,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_a, 	d_a,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_b, 	d_b,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_cxx, 	d_cxx,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_cxy, 	d_cxy,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_cyy, 	d_cyy,	total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_theta, d_theta,total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_abcor, d_abcor,total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_peak,  d_peak, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_dpeak, d_dpeak, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fdpeak,d_fdpeak, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_flux,	d_flux, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_dflux,	d_dflux, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fdflux,	d_fdflux, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fluxerr,	d_fluxerr, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_thresh,	d_thresh, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mthresh,	d_mthresh, total*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fwhm,		d_fwhm, total*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(h_mx,	d_mx, total*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_my,	d_my, total*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mx2,	d_mx2, total*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_my2,	d_my2, total*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mxy,	d_mxy, total*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_poserr_mx2,	d_poserr_mx2, total*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_poserr_my2,	d_poserr_my2, total*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_poserr_mxy,	d_poserr_mxy, total*sizeof(double), cudaMemcpyDeviceToHost);

	cudaMemcpy(h_singuflag,	d_singuflag, total*sizeof(char), cudaMemcpyDeviceToHost);

	for(int i=0; i<NISO; i++)
		cudaMemcpy(h_iso[i],d_iso[i], total*sizeof(int), cudaMemcpyDeviceToHost);

}

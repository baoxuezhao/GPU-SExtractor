/*
 * cudadeblend.h
 *
 *  Created on: 10 Jan, 2013
 *      Author: zhao
 */

#ifndef CUDADEBLEND_H_
#define CUDADEBLEND_H_

void init_deblend(float basethresh,
		int _deblend_nthresh,
		size_t numValidObj,
		size_t numValidPix);

int parcel_out(float basethresh,
		int n_thresh,
		int deb_minarea,
		size_t numValidObj,
		size_t numValidPix,
		double deblend_mincont);

void clear_deblend(int _deblend_nthresh);

#endif /* CUDADEBLEND_H_ */

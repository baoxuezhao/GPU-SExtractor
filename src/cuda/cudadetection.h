/*
 * cudadetection.h
 *
 *  Created on: 10 Jan, 2013
 *      Author: zhao
 */

#ifndef CUDASCAN_H_
#define CUDASCAN_H_

void init_detection();

//void copyScan(float*);

void run_detection(float, float, int, size_t*, size_t*);
void clear_detection();
int compact_and_sort(int numpix, int level, int ext_minarea);

#endif /* CUDASCAN_H_ */

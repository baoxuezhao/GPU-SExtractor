/*
 * cudainit.h
 *
 *  Created on: 10 Jan, 2013
 *      Author: zhao
 */

#ifndef CUDAINIT_H_
#define CUDAINIT_H_

void init_device(int a, int b, float* c);
void clear_device(float *img);

void init_objects(int total);

#endif /* CUDAINIT_H_ */

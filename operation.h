//
//  operation.h
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#ifndef operation_h
#define operation_h

#include <stdio.h>

void vectorAdd(double *v1, double *v2, double *v3);
void vectorSubtract(double *v1, double *v2, double *v3);
void vectorScale(double *v1, const double scale, double *v2);
double vectorNorm(const double *v1);
void vectorNormalize(double *v1, double *v2);
double vectorDotProduct(const double *v1, const double *v2);
void vectorCopy(const double *x1, double *x2);

void directorCalc(double R0, double R, double *position, double *u);

void randinitialize(void);
int randInteger(int max);
double randNormal(void);

#endif /* operation_h */

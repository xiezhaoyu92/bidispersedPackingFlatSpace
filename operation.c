//
//  operation.c
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#include "operation.h"
#include <math.h>
#include "mt19937ar.h"
#include <time.h>


void vectorAdd(double *v1, double *v2, double *v3) {
    v3[0] = v1[0] + v2[0];
    v3[1] = v1[1] + v2[1];
}

void vectorSubtract(double *v1, double *v2, double *v3) {
    v3[0] = v1[0] - v2[0];
    v3[1] = v1[1] - v2[1];
}

void vectorScale(double *v1, const double scale, double *v2) {
    v2[0] = v1[0]*scale;
    v2[1] = v1[1]*scale;
}

double vectorNorm(const double *v1) {
    return sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
}

void vectorNormalize(double *v1, double *v2) {
    double norm = vectorNorm(v1);
    v2[0] = v1[0]/norm;
    v2[1] = v1[1]/norm;
}

double vectorDotProduct(const double *v1, const double *v2) {
    return v1[0]*v2[0] + v1[1]*v2[1];
}

void vectorCopy(const double *x1, double *x2) {
    x2[0] = x1[0];
    x2[1] = x1[1];
}


/***************************
 * random number functions *
 ***************************/



/* initialize random number generator
 * Input:       None
 * Output:      None
 * Returns:     None
 */
void randinitialize(void) {
    FILE *urandom;
    unsigned long seed=(unsigned long) time(NULL);
    char bytes[sizeof(unsigned long)];
    
    /* Initialize the random number generator */
    urandom=fopen("/dev/urandom", "r");
    if (urandom) {
        for(int i=0;i<sizeof(unsigned long);i++) bytes[i]=(char) fgetc(urandom);
        fclose(urandom);
        
        seed = *((unsigned long *) bytes);
    } else printf("Warning: initializing random number generator using time-not recommended for production runs.\n");
    
    init_genrand(seed);
}

/* Generate a random integer from 1 to max
 * Input:	max	      - largest possible number to return
 * Output:      None
 * Returns:	a random integer from 1 to max
 */
int randInteger(int max) {
    
    return (int) (genrand_real1()*max+1.0);
}

/* Generate a random number with a gaussian distribution with sigma=1
 * Input:	None
 * Output:      None
 * Returns:	a random double
 */
double randNormal(void) {
    double x,y,r;
    
    do {
        x=2.0*genrand_real1()-1.0;
        y=2.0*genrand_real1()-1.0;
        
        r=x*x+y*y;
    } while (r>=1.0);
    
    return x*sqrt((-2.0*log(r))/r);
}

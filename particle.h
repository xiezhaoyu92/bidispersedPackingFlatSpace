//
//  particle.h
//  bidispersedPacking
//
//  Created by Zhaoyu Xie on 3/18/19.
//  Copyright © 2019 Zhaoyu Xie. All rights reserved.
//

#ifndef particle_h
#define particle_h

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#include <stdio.h>

typedef struct {
    double position[2];
    double oldPosition[2];
    double postPreRelaxPosition[2];
    
    double force[2];
    double forceStoc[2];
    double rad;
    
    int oldCoord[2];
    int coord[2];
    double postPreRelaxCoord[2];
} particle;

int overlapQ(particle *p1, particle *p2, int nPart[], double L);
int anyOverlapQ(int i, particle p[], int np, int nPart[], double L);
void gradDescentIntegrate(particle *p, double dt, double L);
void partitionOneParticle(particle *p, double partBoundX[], double partBoundY[], int nPart[]);
int addOverlapForce(particle *p1, particle *p2, int nPart[], double L);
void undoOverlaps(particle p[], int maxUndoSteps, double dtOverlap, int *nearJammed, double L, int np, double partBoundX[], double partBoundY[], int nPart[]);
void addStochasticForce(double D, double rad, particle *p);
void updatePartitions(int nPart[], double partBoundX[], double partBoundY[], double L, double rad);
void projectIntoNewArea(particle p[], double L, double LOld);

#endif /* particle_h */

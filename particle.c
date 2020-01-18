//
//  particle.c
//  bidispersedPacking
//
//  Created by Zhaoyu Xie on 3/18/19.
//  Copyright © 2019 Zhaoyu Xie. All rights reserved.
//

#include "particle.h"
#include "operation.h"
#include <math.h>
#include <stdlib.h>

int overlapQ(particle *p1, particle *p2, int nPart[], double L) {
    particle temp;
    temp.rad = p2->rad;
    if(p2->coord[0]==p1->coord[0]||abs(p2->coord[0]-p1->coord[0])==1)
        temp.position[0]=p2->position[0];
    else if(p2->coord[0]-p1->coord[0]==nPart[0]-1)
        temp.position[0]=p2->position[0]-L;
    else if(p2->coord[0]-p1->coord[0]==-(nPart[0]-1))
        temp.position[0]=p2->position[0]+L;
    else
        return FALSE;
    if(p2->coord[1]==p1->coord[1]||abs(p2->coord[1]-p1->coord[1])==1)
        temp.position[1]=p2->position[1];
    else if(p2->coord[1]-p1->coord[1]==nPart[1]-1)
        temp.position[1]=p2->position[1]-L;
    else if(p2->coord[1]-p1->coord[1]==-(nPart[1]-1))
        temp.position[1]=p2->position[1]+L;
    else
        return FALSE;
    
    double r[2];
    vectorSubtract(p1->position, temp.position, r);
    if(vectorNorm(r)-(p1->rad+temp.rad)<MACHINE_EPSILON)
        return TRUE;
    else
        return FALSE;
/*    if((p1->position[0]-temp.position[0])*(p1->position[0]-temp.position[0])+(p1->position[1]-temp.position[1])*(p1->position[1]-temp.position[1])<(p1->rad+temp.rad)*(p1->rad+temp.rad))
        return TRUE;
    else
        return FALSE;*/
}

int anyOverlapQ(int i, particle p[], int np, int nPart[], double L) {
    int overlaps = FALSE;
    for(int j=0; j<np; j++) {
        if(j==i)
            continue;
        if(overlapQ(&(p[i]),&(p[j]),nPart,L)){
            overlaps=TRUE;
            break;
        }
    }
    if(overlaps)
        return TRUE;
    else
        return FALSE;
}

void gradDescentIntegrate(particle *p, double dt, double L) {
    for (int i=0; i<2; i++) {
        p->position[i] = p->position[i] + (p->force[i])*dt + (p->forceStoc[i])*sqrt(dt);
        if(p->position[i]<0)
            p->position[i] = p->position[i]+L;
        if(p->position[i]>L)
            p->position[i] = p->position[i]-L;
    }
}

void partitionOneParticle(particle *p, double partBoundX[], double partBoundY[], int nPart[]) {
    int i=0, j=0;
    while(i<nPart[0]&&p->position[0]>=partBoundX[i])
        i++;
    p->coord[0] = i-1;
    while(j<nPart[1]&&p->position[1]>=partBoundY[j])
        j++;
    p->coord[1] = j-1;
}

int addOverlapForce(particle *p1, particle *p2, int nPart[], double L) {
    if(overlapQ(p1,p2,nPart,L)) {
        particle temp;
        temp.rad = p2->rad;
        temp.force[0] = 0;
        temp.force[1] = 0;
        if(p2->coord[0]==p1->coord[0]||abs(p2->coord[0]-p1->coord[0])==1)
            temp.position[0]=p2->position[0];
        else if(p2->coord[0]-p1->coord[0]==nPart[0]-1)
            temp.position[0]=p2->position[0]-L;
        else if(p2->coord[0]-p1->coord[0]==-(nPart[0]-1))
            temp.position[0]=p2->position[0]+L;
        if(p2->coord[1]==p1->coord[1]||abs(p2->coord[1]-p1->coord[1])==1)
            temp.position[1]=p2->position[1];
        else if(p2->coord[1]-p1->coord[1]==nPart[1]-1)
            temp.position[1]=p2->position[1]-L;
        else if(p2->coord[1]-p1->coord[1]==-(nPart[1]-1))
            temp.position[1]=p2->position[1]+L;
        double deltaX[2];
        double distance;
        vectorSubtract(p1->position, temp.position, deltaX);
        distance = vectorNorm(deltaX);
        for(int i=0; i<2; i++){
            p1->force[i] += (p1->rad+temp.rad)*deltaX[i]/distance;
            temp.force[i] = -(p1->rad+temp.rad)*deltaX[i]/distance;
            p2->force[i] += temp.force[i];
        }
        return TRUE;
    }
    else
        return FALSE;
}

void undoOverlaps(particle p[], int maxUndoSteps, double dtOverlap, int *nearJammed, double L, int np, double partBoundX[], double partBoundY[], int nPart[]) {
    for(int i=0; i<np; i++)
        for(int j=0; j<2; j++)
            p[i].forceStoc[j] = 0;
    
    int undoSteps = 0;
    *nearJammed = FALSE;
    int totalOverlapQ;
    do {
        totalOverlapQ = FALSE;
        undoSteps++;
        //printf("%d\n", undoSteps);
        if(undoSteps>maxUndoSteps) {
            *nearJammed = TRUE;
            break;
        }
        for(int i=0; i<np; i++)
            for(int j=0; j<2; j++)
                p[i].force[j] = 0;
        for(int i=0; i<np; i++)
            for(int j=i+1; j<np; j++){
                if(addOverlapForce(&(p[i]),&(p[j]),nPart,L))
                    totalOverlapQ = TRUE;
            }
        for(int i=0; i<np; i++) {
            vectorCopy(p[i].position, p[i].oldPosition);
            gradDescentIntegrate(&(p[i]), dtOverlap, L);
            if(isnan(p[i].position[0]))
                vectorCopy(p[i].oldPosition, p[i].position);
            partitionOneParticle(&(p[i]), partBoundX, partBoundY, nPart);
        }
    } while(totalOverlapQ);
}

//rad is the radius of larger sphere
void addStochasticForce(double D, double rad, particle *p) {
    for(int i=0; i<2; i++)
        p->forceStoc[i] = sqrt(2*D)/rad*(p->rad)*randNormal();
}

//calculate the total number of partitions
void updatePartitions(int nPart[], double partBoundX[], double partBoundY[], double L, double rad){
    nPart[0] = (int)(L/(2*rad));
    nPart[1] = (int)(L/(2*rad));
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = i*2*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = i*2*rad;
}

void projectIntoNewArea(particle *p, double L, double LOld) {
        for(int j=0; j<2; j++)
            p->position[j] = p->position[j]/LOld*L;
}

//
//  main.c
//  bidispersedPacking
//
//  Created by Zhaoyu Xie on 3/15/19.
//  Copyright © 2019 Zhaoyu Xie. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <stdlib.h>
#include "particle.h"
#include "operation.h"
#include "mt19937ar.h"

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON = 1e-15
#endif

int main(int argc, char * argv[]) {
    double L0 = 1; //initial size of area
    int np = 100; //number of particles
    double rad = 0.01; //raidus of large disks
    double mixRatio = 0.5; //ratio of number of large disks to small ones
    double ratio = 1; //ratio of small radius to large radius
    
    double diffusionTimescale = 0.25;
    
    int animationQ = FALSE;
    int animationRate = 100;
    int readUnfinishedQ = FALSE;
    int readJammedQ = FALSE;
    
    FILE* acceptFile = NULL; //record the acceptance rate and dtDiffusion
    int acceptRecordQ = TRUE;
    FILE *animationFile=NULL;
    
    static struct option longOptions[]={
        {"radius",             required_argument,  NULL, 'r'},
        {"radRatio",           required_argument,  NULL, 'R'},
        {"particleNumber",     required_argument,  NULL, 'n'},
        {"mixRatio",           required_argument,  NULL, 'i'},
        {"initialSize",        required_argument,  NULL, 'l'},
        {"animate",            no_argument,        NULL, 'm'},
        {"dtDiffusionScale",   required_argument,  NULL, 's'},
        {"readUnfinished",     no_argument,        NULL, 'u'},
        {"readJammed",         no_argument,        NULL, 'j'},
        {0,                    0,                  0,     0 }
    };
    int optIndex = 0;
    int opt;
    while ((opt = getopt_long(argc, argv,"r:R:n:i:l:ms:uj",
                              longOptions, &optIndex )) != -1) {
        switch (opt) {
            case 0:
                break;
            case 'r' :
                rad = atof(optarg);
                break;
            case 'R' :
                ratio = atof(optarg);
                break;
            case 'n' :
                np = atoi(optarg);
                break;
            case 'i' :
                mixRatio = atof(optarg);
                break;
            case 'l' :
                L0 = atof(optarg);
                break;
            case 'm' :
                animationQ = TRUE;
                break;
            case 's':
                diffusionTimescale = atof(optarg);
                break;
            case 'u' :
                readUnfinishedQ = TRUE;
                break;
            case 'j' :
                readJammedQ = TRUE;
                break;
            default:
                exit(1);
        }
    }
    
    if (readUnfinishedQ && readJammedQ) {
        printf("readJammmed and readUnfinished are mutually exclusive options. Please specify only one.\n");
        exit(1);
    }
    
    if(readUnfinishedQ || readJammedQ) {
        FILE *npFile=NULL;
        npFile = fopen("npts.dat","r");
        if (npFile) {
            fscanf(npFile, "%i", &np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
            exit(1);
        }
        FILE *radFile=NULL;
        radFile = fopen("rad.dat","r");
        if (radFile) {
            fscanf(radFile, "%lf", &rad);
            fclose(radFile);
        } else {
            printf("radFile pointer is null\n");
            exit(1);
        }
    }
    else {
        FILE *npFile=NULL;
        npFile = fopen("npts.dat","w");
        if (npFile) {
            fprintf(npFile, "%i\n", np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
            exit(1);
        }
        FILE *radFile=NULL;
        radFile = fopen("rad.dat","w");
        if (radFile) {
            fprintf(radFile, "%lf\n", rad);
            fclose(radFile);
        } else {
            printf("radFile pointer is null\n");
            exit(1);
        }
        FILE *ratioFile=NULL;
        ratioFile = fopen("radRatio.dat","w");
        if (ratioFile) {
            fprintf(ratioFile, "%lf\n", ratio);
            fclose(ratioFile);
        } else {
            printf("ratioFile pointer is null\n");
            exit(1);
        }
        FILE *mixRatioFile=NULL;
        mixRatioFile = fopen("mixRatio.dat","w");
        if (mixRatioFile) {
            fprintf(mixRatioFile, "%lf\n", mixRatio);
            fclose(mixRatioFile);
        } else {
            printf("mixRatioFile pointer is null\n");
            exit(1);
        }
        FILE *initSizeFile=NULL;
        initSizeFile = fopen("initialSize.dat","w");
        if (initSizeFile) {
            fprintf(initSizeFile, "%lf\n", L0);
            fclose(initSizeFile);
        } else {
            printf("initSizeFile pointer is null\n");
            exit(1);
        }
    }
    
    particle p[np];
    
    double L;
    double LOld;
    
    int jammed = FALSE;
    int maxUndoStepsStart = 1e6;
    double dtOverlap = 1e-4;
    double dtOverlapStartup = 1e-4;
    int maxUndoSteps = 1e4;
    
    double simTimeStart = 0;
    double nextRelaxationTime = 0;
    int relaxationStep = 0;
    int simStep = 0;
    double simTime;
    
    double diffusionCoeff = pow(2*rad/1000,2)/2;
    double simTimeFinal = 2*rad*rad/(diffusionCoeff*diffusionTimescale);
    double relaxationSteps = 1e5;
    double dtRelaxation = simTimeFinal/relaxationSteps;
    double dtTol = 1e-5*dtRelaxation;
    int outputSteps = 10000;
    
    double dtDiffusion = 1;
    //double dtDiffusionMin = 1e-3;
    double timestepScaleFactor = 100; // how much to scale the preliminary timestep, which gives an average step equal to the interparticle spacing and results in an acceptance ratio of ~0.5
    double minSpacing = 1e-7; // spacing on which to base minimum dtDiffusion
    double dtDiffusionMin = timestepScaleFactor*minSpacing*minSpacing/(2*diffusionCoeff);
    if (dtDiffusion>dtRelaxation) {
        printf("Fast relaxation rate: using smaller timestep.\n");
        dtDiffusion = dtRelaxation;
        printf("%e\n",dtDiffusion);
    }
    
    if(readUnfinishedQ||readJammedQ) {
        FILE *initSizeFile=NULL;
        initSizeFile = fopen("initialSize.dat","r");
        if (initSizeFile) {
            fscanf(initSizeFile, "%lf", &L0);
            fclose(initSizeFile);
        } else {
            printf("initSizeFile pointer is null\n");
            exit(1);
        }
        FILE *radiusFile=NULL;
        radiusFile = fopen("radius.dat","r");
        if (radiusFile) {
            for(int i=0; i<np; i++)
                fscanf(radiusFile, "%lf", &(p[i].rad));
            fclose(radiusFile);
        } else {
            printf("radiusFile pointer is null\n");
            exit(1);
        }
        if(readUnfinishedQ){
            FILE *sizeTemp=NULL;
            sizeTemp = fopen("sizeTemp.dat","r");
            if (sizeTemp) {
                fscanf(sizeTemp, "%lf", &L);
                fclose(sizeTemp);
            } else {
                printf("sizeTemp pointer is null\n");
                exit(1);
            }
            FILE *simTimeTempFile=NULL;
            simTimeTempFile = fopen("simTimeTemp.dat","r");
            if (simTimeTempFile) {
                fscanf(simTimeTempFile, "%lf", &simTimeStart);
                fclose(simTimeTempFile);
            } else {
                printf("simTimeTempFile pointer is null\n");
                exit(1);
            }
            FILE *nextRelaxationFile=NULL;
            nextRelaxationFile = fopen("nextRelaxationTime.dat","r");
            if (nextRelaxationFile) {
                fscanf(nextRelaxationFile, "%lf", &nextRelaxationTime);
                fclose(nextRelaxationFile);
            } else {
                printf("nextRelaxationFile pointer is null\n");
                exit(1);
            }
            FILE *steps = NULL;
            steps = fopen("stepsTemp.dat","r");
            if(steps) {
                fscanf(steps, "%d %d", &simStep, &relaxationStep);
                fclose(steps);
            }
            else {
                printf("steps pointer is null\n");
                exit(1);
            }
            FILE *dtFile=NULL;
            dtFile = fopen("dtTemp.dat","r");
            if (dtFile) {
                fscanf(dtFile, "%lf %lf %lf", &dtDiffusion, &dtRelaxation, &dtOverlap);
                fclose(dtFile);
            } else {
                printf("dtFile pointer is null\n");
                exit(1);
            }
        }
        else {
            FILE *size=NULL;
            size = fopen("size.dat","r");
            if (size) {
                fscanf(size, "%lf", &L);
                fclose(size);
            } else {
                printf("size pointer is null\n");
                exit(1);
            }
            FILE *simTimeFile=NULL;
            simTimeFile = fopen("simTimeArrested.dat","r");
            if (simTimeFile) {
                fscanf(simTimeFile, "%lf", &simTimeStart);
                fclose(simTimeFile);
            } else {
                printf("simTimeTempFile pointer is null\n");
                exit(1);
            }
            nextRelaxationTime = simTimeStart;
        }
    }
    else {
        L = L0;
        for (int i=0; i< (int)(np*mixRatio); i++) {
            p[i].rad=rad;
        }
        for (int i=(int)(np*mixRatio); i<np; i++) {
            p[i].rad=rad*ratio;
        }
        FILE *radiusFile=NULL;
        radiusFile = fopen("radius.dat","w");
        if (radiusFile) {
            for(int i=0; i<np; i++)
                fprintf(radiusFile, "%lf\n", p[i].rad);
            fclose(radiusFile);
        } else {
            printf("radiusFile pointer is null\n");
            exit(1);
        }
    }
    
    double simTimeLastRelax = simTimeStart;
    double nextRelaxationTimeLastRelax = nextRelaxationTime;
    
    int nPart[2];
    nPart[0] = (int)(L/(2*rad));
    nPart[1] = (int)(L/(2*rad));
    double partBoundX[nPart[0]];
    double partBoundY[nPart[1]];
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = i*2*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = i*2*rad;
    
    if(readUnfinishedQ) {
        FILE *configTemp=NULL;
        configTemp = fopen("configurationTemp.asc","r");
        if (configTemp) {
            for(int i=0; i<np; i++){
                fscanf(configTemp, "%lf %lf", &(p[i].position[0]), &(p[i].position[1]));
                partitionOneParticle(&(p[i]), partBoundX, partBoundY, nPart);
                vectorCopy(p[i].position, p[i].postPreRelaxPosition);
                for (int l=0; l<2; l++)
                    p[i].postPreRelaxCoord[l] = p[i].coord[l];
            }
            fclose(configTemp);
        } else {
            printf("configTemp pointer is null\n");
            exit(1);
        }
    }
    else if(readJammedQ){
        FILE *config=NULL;
        config = fopen("configuration.asc","r");
        if (config) {
            for(int i=0; i<np; i++) {
                fscanf(config, "%lf %lf", &(p[i].position[0]), &(p[i].position[1]));
                partitionOneParticle(&(p[i]), partBoundX, partBoundY, nPart);
            }
            fclose(config);
        } else {
            printf("config pointer is null\n");
            exit(1);
        }
        double dtOverlapMin = 1e-5;
        undoOverlaps(p, maxUndoStepsStart, dtOverlapMin, &jammed, L, np, partBoundX, partBoundY, nPart);
        if (jammed) {
            printf("overcrowded starting point - there is overlap in the last configuration\n");
            exit(1);
        }
        for(int j=0; j<np; j++) {
            vectorCopy(p[j].position, p[j].postPreRelaxPosition);
            for (int l=0; l<2; l++)
                p[j].postPreRelaxCoord[l] = p[j].coord[l];
        }
    }
    else {
        
        randinitialize();
        
        FILE *init1=NULL;
        init1 = fopen("init1.asc","w");
        if(!init1) {
            printf("file init1.asc failed to open");
            exit(1);
        }
        for(int i=0; i<np; i++) {
            double x = L*genrand_real2();
            double y = L*genrand_real2();
            p[i].position[0] = x;
            p[i].position[1] = y;
            fprintf(init1, "%.15lf %.15lf\n", x, y);
            partitionOneParticle(&(p[i]), partBoundX, partBoundY, nPart);
        }
        fclose(init1);
        
        undoOverlaps(p, maxUndoStepsStart, dtOverlap, &jammed, L, np, partBoundX, partBoundY, nPart);
        if (jammed) {
            printf("overcrowded starting point - use fewer particles or smaller radius\n");
            exit(1);
        }
        FILE *init2=NULL; // after relaxation
        init2 = fopen("init2.asc","w");
        if (!init2) {
            printf("file init2.asc failed to open");
            exit(1);
        }
        for (int i = 0; i < np; i++) {
            fprintf(init2, "%.15lf %.15lf\n", p[i].position[0], p[i].position[1]);
        }
        fclose(init2);
        
        FILE *init3=NULL;;
        init3 = fopen("init3.asc","w");
        if (!init3) {
            printf("file init3.asc failed to open");
            exit(1);
        }
        for(int i=0; i<np; i++)
            for(int j=0; j<2; j++)
                p[i].force[j]=0;
        int startupDiffusionSteps = 20000;
        double diffusionCoeffStartup = pow(2*rad/100, 2);
        for (int i=0; i<startupDiffusionSteps; i++) {
            if(i%1000==0) printf("startup shuffle: %i\n",i);
            for (int j=0; j<np; j++) {
                int k = randInteger(np)-1;
                vectorCopy(p[k].position, p[k].oldPosition);
                for (int l=0; l<2; l++)
                    p[k].oldCoord[l] = p[k].coord[l];
                addStochasticForce(diffusionCoeffStartup, rad, &(p[k]));
                gradDescentIntegrate(&(p[k]), dtDiffusion, L);
                partitionOneParticle(&(p[k]), partBoundX, partBoundY, nPart);
                if(anyOverlapQ(k, p, np, nPart, L)||isnan(p[k].position[0])) {
                    vectorCopy(p[k].oldPosition, p[k].position);
                    for (int l=0; l<2; l++)
                        p[k].coord[l] = p[k].oldCoord[l];
                }
            }
            /*if (i%animationRate==0)
                for (int j = 0; j < np; j++) {
                    fprintf(init3, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
                }*/
        }
        for (int j = 0; j < np; j++) {
            fprintf(init3, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
        }
        fclose(init3);
        
        if (animationQ) {
            animationFile = fopen("animation.dat", "a");
            if (animationFile) {
                fprintf(animationFile, "%.15lf\n", L);
                for (int j=0; j<np; j++)
                    fprintf(animationFile, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
                fclose(animationFile);
            }
            else
                printf("animationFile pointer is null\n");
        }

    }
    
    
    for (simTime = simTimeStart; simTime <= simTimeFinal && !jammed; simTime += dtDiffusion) {
        if (simTime > nextRelaxationTime) {
            LOld = L;
            for(int j=0; j<np; j++) {
                vectorCopy(p[j].position, p[j].oldPosition);
                for (int l=0; l<2; l++)
                    p[j].oldCoord[l] = p[j].coord[l];
            }
            L = L0*(1-simTime/simTimeFinal);
            updatePartitions(nPart, partBoundX, partBoundY, L, rad);
            for(int j=0; j<np; j++) {
                projectIntoNewArea(&(p[j]), L, LOld);
                partitionOneParticle(&(p[j]), partBoundX, partBoundY, nPart);
            }
            int nearJammed = FALSE;
            undoOverlaps(p, maxUndoSteps, dtOverlap, &nearJammed, L, np, partBoundX, partBoundY, nPart);
            if(nearJammed) {
                L = LOld;
                for(int j=0; j<np; j++) {
                    vectorCopy(p[j].postPreRelaxPosition, p[j].position);
                    for (int l=0; l<2; l++)
                        p[j].coord[l] = p[j].postPreRelaxCoord[l];
                }
                updatePartitions(nPart, partBoundX, partBoundY, L, rad);
                simTime = simTimeLastRelax;
                nextRelaxationTime = nextRelaxationTimeLastRelax;
                dtRelaxation *= 0.5;
                dtOverlap *= 0.5;
                printf("reducing timesteps %e %e\n",dtRelaxation,dtOverlap);
                if (dtRelaxation < dtTol)
                    jammed = TRUE;
                if (dtDiffusion > dtRelaxation)
                    dtDiffusion = dtRelaxation;
            }
            else {
                simTimeLastRelax = simTime;
                nextRelaxationTimeLastRelax = nextRelaxationTime;
                for(int j=0; j<np; j++) {
                    vectorCopy(p[j].position, p[j].postPreRelaxPosition);
                    for (int l=0; l<2; l++)
                        p[j].postPreRelaxCoord[l] = p[j].coord[l];
                }
                relaxationStep++;
                if (animationQ && relaxationStep % animationRate == 0) {
                    animationFile = fopen("animation.dat", "a");
                    if (animationFile) {
                        fprintf(animationFile, "%.15lf\n", L);
                        for (int j=0; j<np; j++)
                            fprintf(animationFile, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
                        fclose(animationFile);
                    } else
                        printf("animationFile pointer is null\n");
                }
            }
            nextRelaxationTime += dtRelaxation;
        }
        
        //***************************//
        // one set of particle moves //
        //***************************//
        int stepCount = 0;
        int overlapCount = 0;
        for(int i=0; i<np; i++)
            for(int j=0; j<2; j++){
                p[i].force[j] = 0;
                p[i].forceStoc[j] = 0;
            }
        
        for(int i=0; i<np; i++)
            addStochasticForce(diffusionCoeff, rad, &(p[i]));
        // going to loop through particles in a random order - set up shuffled array of indices
        int randIndices[np];
        for (int k=0; k<np; k++) {
            randIndices[k]=k;
        }
        for (int k=np-1; k>0; k--) {
            int swapWith = randInteger(k-1);
            int hold = randIndices[k];
            randIndices[k] = randIndices[swapWith];
            randIndices[swapWith] = hold;
        }
        
        for (int k0 = 0; k0 < np; k0++) {
            int k = randIndices[k0];
            vectorCopy(p[k].position, p[k].oldPosition);
            for (int l=0; l<2; l++)
                p[k].oldCoord[l] = p[k].coord[l];
            gradDescentIntegrate(&(p[k]), dtDiffusion, L);
            partitionOneParticle(&(p[k]), partBoundX, partBoundY, nPart);
            stepCount++;
            if(anyOverlapQ(k, p, np, nPart, L)||isnan(p[k].position[0])) {
                vectorCopy(p[k].oldPosition, p[k].position);
                for (int l=0; l<2; l++)
                    p[k].coord[l] = p[k].oldCoord[l];
                overlapCount++;
            }
        }
        double acceptanceRatio = 1.0 - (double)overlapCount/stepCount;
        
        if (acceptRecordQ) {
            acceptFile = fopen("acceptance.dat","a");
            if (acceptFile) {
                fprintf(acceptFile, "%e %f %e\n", dtDiffusion, acceptanceRatio, dtRelaxation);
                fclose(acceptFile);
            } else {
                printf("acceptFile pointer is null\n");
            }
        }
        
        if(acceptanceRatio < 0.5)
            dtDiffusion = dtDiffusion*0.99;
        else
            dtDiffusion = dtDiffusion*1.01;
        if(dtDiffusion < dtDiffusionMin)
            dtDiffusion = dtDiffusionMin;
        if(dtDiffusion > dtRelaxation)
            dtDiffusion = dtRelaxation;
        
        if (simStep%outputSteps==0) {
            FILE *sizeTemp=NULL;
            sizeTemp = fopen("sizeTemp.dat","w");
            if (sizeTemp) {
                fprintf(sizeTemp, "%.15lf\n", L);
                fclose(sizeTemp);
            } else {
                printf("sizeTemp pointer is null\n");
                exit(1);
            }
            FILE *configTemp=NULL;
            configTemp = fopen("configurationTemp.asc","w");
            if (configTemp) {
                for(int i=0; i<np; i++)
                    fprintf(configTemp, "%.15lf %.15lf\n", p[i].position[0], p[i].position[1]);
                fclose(configTemp);
            } else {
                printf("configTemp pointer is null\n");
                exit(1);
            }
            FILE *simTimeTempFile=NULL;
            simTimeTempFile = fopen("simTimeTemp.dat","w");
            if (simTimeTempFile) {
                fprintf(simTimeTempFile, "%.15lf\n", simTime);
                fclose(simTimeTempFile);
            } else {
                printf("simTimeTempFile pointer is null\n");
                exit(1);
            }
            FILE *nextRelaxationFile=NULL;
            nextRelaxationFile = fopen("nextRelaxationTime.dat","w");
            if (nextRelaxationFile) {
                fprintf(nextRelaxationFile, "%.15lf\n", nextRelaxationTime);
                fclose(nextRelaxationFile);
            } else {
                printf("nextRelaxationFile pointer is null\n");
                exit(1);
            }
            FILE *steps = NULL;
            steps = fopen("stepsTemp.dat","w");
            if(steps) {
                fprintf(steps, "%d %d\n", simStep+1, relaxationStep);
                fclose(steps);
            }
            else {
                printf("steps pointer is null\n");
                exit(1);
            }
            FILE *dtFile=NULL;
            dtFile = fopen("dtTemp.dat","w");
            if (dtFile) {
                fprintf(dtFile, "%.15lf %.15lf %.15lf\n", dtDiffusion, dtRelaxation, dtOverlap);
                fclose(dtFile);
            } else {
                printf("dtFile pointer is null\n");
                exit(1);
            }
        }
        simStep++;
        printf("simTime/simTimeFinal: %.16lf\n", simTime/simTimeFinal);
    }
    
    FILE *sizeFile=NULL;
    sizeFile = fopen("size.dat","w");
    if (sizeFile) {
        fprintf(sizeFile, "%.15lf\n", L);
        fclose(sizeFile);
    } else {
        printf("sizeFile pointer is null\n");
        exit(1);
    }
    FILE *configFile=NULL;
    configFile = fopen("configuration.asc","w");
    if (configFile) {
        for(int i=0; i<np; i++)
            fprintf(configFile, "%.15lf %.15lf\n", p[i].position[0], p[i].position[1]);
        fclose(configFile);
    } else {
        printf("configFile pointer is null\n");
        exit(1);
    }
    FILE *simTimeArrestedFile=NULL;
    simTimeArrestedFile = fopen("simTimeArrested.dat","w");
    if (simTimeArrestedFile) {
        fprintf(simTimeArrestedFile, "%.15lf\n", simTime);
        fclose(simTimeArrestedFile);
    } else {
        printf("simTimeArrestedFile pointer is null\n");
        exit(1);
    }
    FILE *nextRelaxationArrestedFile=NULL;
    nextRelaxationArrestedFile = fopen("nextRelaxationTimeArrested.dat","w");
    if (nextRelaxationArrestedFile) {
        fprintf(nextRelaxationArrestedFile, "%.15lf\n", nextRelaxationTime);
        fclose(nextRelaxationArrestedFile);
    } else {
        printf("nextRelaxationArrestedFile pointer is null\n");
        exit(1);
    }
    FILE *stepsFile = NULL;
    stepsFile = fopen("steps.dat","w");
    if(stepsFile) {
        fprintf(stepsFile, "%d %d\n", simStep, relaxationStep);
        fclose(stepsFile);
    }
    else {
        printf("stepsFile pointer is null\n");
        exit(1);
    }

    printf("Hello, World!\n");
    return 0;
    
}

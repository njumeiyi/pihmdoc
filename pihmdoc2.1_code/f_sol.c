#define DEBUG
/*******************************************************************************
 *-----------------------------------------------------------------------------*
 * File        : f.c   (PIHM v.2.0)                                            *
 * Function    : Model Kernel: Building ODE system for each physical process   *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * Developer of PIHM v.2.0:  Mukesh Kumar (muk139@psu.edu)		       *
 * Developer of PIHM v.1.0:  Yizhong Qu	(quyizhong@gmail.com)		       *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * NOTE: f.c has gone a massive revamp (essentially rewritten) since PIHM v.1.0*
 *                                                                             *
 *...........MODFICATIONS/ADDITIONS incorporated in f.c (PIHM v.2.0)...........*
 * a) Surface Flow: 							       *
 *	--> Correction of diffusion wave approximation (calculation of dh/ds)  *
 *              i. Calculation of dh/ds performed using planar slope connecting*
 *                 neighboring centroids				       *
 *              ii.Reflection of elements at boundaries and rivers for dh/ds   *
 *		   calculation
 *	--> Correction of kinematic wave approximation (dh/ds calculation based*
 *	    on elevation only instead of head				       *
 *	--> Correction of gradient for cases with steep change in topography   *
 * b) Subsurface Flow:							       *
 *	--> Addition of macropore phenomena				       *
 *	--> Addition of rectangular cell beneath a river element	       *
 *	--> Implementation of two layered subsurface model(sat/unsat) based on *
 *	Richard's eqn							       *
 *	--> Incorporation of Vertical and Horizontal Anisotropy                *
 *	--> Use of geologic data					       *
 * c) River Flow:							       *
 *	--> Correction of kinematic and diff. wave approximation of SV eqn     *
 *	--> Incorporation of flexible river shapes			       *
 *	--> Separate incorporation of leakage and lateral flow		       *
 *	--> Correction of bank overland flow for extreme cases		       *
 *	--> Addition of aquifer cells below river elements		       *
 * c) Surface/Subsurface Coupling:					       *
 *	--> Implementation of First order coupling through (in/ex)filtration   *
 *		based on head continuity 				       *
 * d) Evaporation:							       *
 *	--> Incorporation of ET from ground/subsurface/vegetation	       *
 *	--> Incorporation of landcover properties for calculation of each ET   *
 *	    component							       *
 * e) Computational:							       *
 *	--> Use of temporary state variables in calculation. Note: Never change* 
 *		core state variables					       *
 * f) Miscellaneous (other advantages realtive to PIHM1.0): No maximum         *
 *    constraint on gw level. Accordingly, no numerical constraints on subsur- *
 *    face flux terms.Faster Implementation. Led to first large scale model    *
 *    application.
 *-----------------------------------------------------------------------------*
 *									       *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact 				       *
 *	--> Mukesh Kumar (muk139@psu.edu)				       *
 *	--> Prof. Chris Duffy (cxd11@psu.edu)				       *
 * This code is free for research purpose only.		       		       *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * REFERENCES:                                                                 *
 * PIHM2.0:                                                                    *
 *      a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *      Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *      b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *      a Mesoscale Watershed", Advances in Water Resources (submitted)        *
 * PIHM1.0:                                                                    *
 *      a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *      simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *      b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *      for multiprocess watershed simulation". Water Resources Research       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundials_types.h"   
#include "pihm.h"

#include "f.h"
      
#define multF	2
#define MINpsi	-70
#define EPS 0.05
#define THRESH 0.0
#define UNIT_C 1440	/* Note 60*24 for calculation of yDot in m/min units while forcing is in m/day. */
#define GRAV 9.8*60*60	/* Note the dependence on physical units */ 

#define C_air 1004.0
#define Lv (2.503*pow(10,6))
#define SIGMA (5.67*pow(10,-8)*60)
#define R_dry 287.04
#define R_v 461.5

/* TRANSPORT DEFINES - START */
#define AGE if(MD->AGE == 1)
#define PI 3.14159265
#define ABS_TOL EPS/100
/* TRANSPORT DEFINES - END */



/* TRANSPORT: Initialize DC terms - START */
void initDummyDC(realtype *DummyDC, Model_Data MD){
    int i, j;
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<3*MD->NumEle; j++){
             DummyDC[i*3*MD->NumEle + j] = 0.0;
        }
    }
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<2*MD->NumRiv; j++){
             DummyDC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = 0.0;
        }
    }
	/* Two test-numerical experiments */
	/*
	for(i=0; i<1; i++){
		for(j=1*MD->NumEle; j<2*MD->NumEle; j++){
			DummyDC[i*3*MD->NumEle + j] = 0.000000001;
		}
	}
	*/
	/*
	for(i=0; i<1; i++){
		for(j=0; j<MD->NumRiv; j++){
			DummyDC[(MD->NumSolute*3*MD->NumEle) + i*MD->NumRiv + j] = 0.1;
		}
	}
	*/
}
/* TRANSPORT: Initialize DC terms - END */

/* TRANSPORT: Initialize C terms - START */
void initDummyC(realtype *DummyC, Model_Data MD, realtype *Y, realtype t, realtype C){
    int i, j;
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<3*MD->NumEle; j++){
//??		Y[(3*MD->NumEle+MD->NumRiv) + i*3*MD->NumEle + j] = Y[j]>(j<2*MD->NumEle?ABS_TOL:0) ? Y[(3*MD->NumEle+MD->NumRiv) + i*3*MD->NumEle + j] : 0.0;
	     //??if(Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] > 0 && Y[j] > EPS/100)
	     //TODO if(Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] > 0 && Y[j] > 0)
		 //if(Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] > 0)
		 if(MD->EleConc[j][i]>0)
		     {DummyC[i*3*MD->NumEle + j] = Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] / (Y[j]*MD->Ele[j%MD->NumEle].area); //Mass to concentration (Mass/Volume)	
		     //DummyC[i*3*MD->NumEle + j] = Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j];//?? TODO
			DummyC[i*3*MD->NumEle + j] = MD->EleConc[j][i];
	//if(j==2*MD->NumEle+9) {printf("\n->%lf\n", DummyC[i*3*MD->NumEle + j]); getchar();}
 }
	             //OLD DummyC[i*3*MD->NumEle + j] = Y[j]>ABS_TOL?Y[(3*MD->NumEle+MD->NumRiv) + i*3*MD->NumEle + j]:0.0;
	     else
		     //??DummyC[i*3*MD->NumEle + j] = 0.0;
			DummyC[i*3*MD->NumEle + j] = C; //TODO FIX ME: Interpolation(&MD->TSD_PChem[i*MD->NumPrep+MD->Ele[j].prep-1], t);
        }
    }
//printf("\n@");
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<2*MD->NumRiv; j++){
//??		Y[(3*MD->NumEle+MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j] = Y[3*MD->NumEle+j]>0 ? Y[(3*MD->NumEle+MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j] : 0.0;
	     if(Y[(3*MD->NumEle+2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] > 0)
	            //OLD DummyC[(MD->NumSolute*3*MD->NumEle) + i*MD->NumRiv + j] = Y[3*MD->NumEle+j]>ABS_TOL?Y[(3*MD->NumEle+MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j]:0.0;
	//??	    DummyC[(MD->NumSolute*3*MD->NumEle) + i*MD->NumRiv + j] = Y[(3*MD->NumEle+2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j] / (Y[3*MD->NumEle+j]*MD->Riv[j%MD->NumRiv].coeff*MD->Riv[j%MD->NumRiv].Length);
		     //DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = Y[(3*MD->NumEle+2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]; //?? TODO
		    DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = MD->RivConc[j][i];
	     else
		     DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = 0.0;
		if(j==0*MD->NumRiv + 5) printf(" %lf\t",DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j]);
        }
    }
	//??
	/*
    if(t<=10){
    	DummyC[0 + 2*MD->NumEle+667-1] = 1.0; //??
	DummyC[667-1] = 2.0; //??
	DummyC[MD->NumSolute*3*MD->NumEle+0 + 0*MD->NumRiv] = 3.0; //??
    }
    else{
	DummyC[0 + 2*MD->NumEle+667-1] = 0.0; //??
        DummyC[667-1] = 0.0; //??
        DummyC[MD->NumSolute*3*MD->NumEle+0 + 0*MD->NumRiv] = 0.0; //??
    }
	//??*/
}
/* TRANSPORT: Initialize C terms - END */


realtype SLOPE(realtype x1, realtype y1, realtype x2, realtype y2){
        double delX,delY;
        delX = x2 - x1;
        delY = (delX==0?y2 - y1 + 0.0000001:y2 - y1);
        if(delX<0.0 && delY>0.0)
                return PI - atan(fabs(delY)/fabs(delX));
        else if(delX<=0.0 && delY<=0.0)
                return PI + atan(fabs(delY)/fabs(delX));
        else if(delX>0.0 && delY<0.0)
                return 2*PI - atan(fabs(delY)/fabs(delX));
        else
                return 0 + atan(fabs(delY)/fabs(delX));
}


/* AGE: Initialize AgeDC terms - END */
void initDummyAgeDC(realtype *DummyAgeDC, Model_Data MD){
    int i, j;
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<3*MD->NumEle; j++){
             DummyAgeDC[i*3*MD->NumEle + j] = 0.0;
        }
    }
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<2*MD->NumRiv; j++){
             DummyAgeDC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = 0.0;
        }
    }
}

/* AGE: Initialize AgeDC terms - END */

/* AGE: Initialize AgeC terms - END */
void initDummyAgeC(realtype *DummyAgeC, Model_Data MD, realtype *Y, realtype t){

    int i, j;
	realtype C; C = 0.0;
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<3*MD->NumEle; j++){
		 if(MD->EleAgeConc[j][i]>0)
		     {DummyAgeC[i*3*MD->NumEle + j] = Y[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j] / (Y[j]*MD->Ele[j%MD->NumEle].area); //Mass to concentration (Mass/Volume)	
		     //DummyC[i*3*MD->NumEle + j] = Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j];//?? TODO
			DummyAgeC[i*3*MD->NumEle + j] = MD->EleAgeConc[j][i];
	//if(j==2*MD->NumEle+9) {printf("\n->%lf\n", DummyC[i*3*MD->NumEle + j]); getchar();}
 }
	             //OLD DummyC[i*3*MD->NumEle + j] = Y[j]>ABS_TOL?Y[(3*MD->NumEle+MD->NumRiv) + i*3*MD->NumEle + j]:0.0;
	     else
		     //??DummyC[i*3*MD->NumEle + j] = 0.0;
			DummyAgeC[i*3*MD->NumEle + j] = C; //TODO FIX ME: Interpolation(&MD->TSD_PChem[i*MD->NumPrep+MD->Ele[j].prep-1], t);
        }
    }
//printf("\n@");
    for(i=0; i<MD->NumSolute; i++){
        for(j=0; j<2*MD->NumRiv; j++){
//??		Y[(3*MD->NumEle+MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j] = Y[3*MD->NumEle+j]>0 ? Y[(3*MD->NumEle+MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j] : 0.0;
	     if(Y[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] > 0)
	            //OLD DummyC[(MD->NumSolute*3*MD->NumEle) + i*MD->NumRiv + j] = Y[3*MD->NumEle+j]>ABS_TOL?Y[(3*MD->NumEle+MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j]:0.0;
	//??	    DummyC[(MD->NumSolute*3*MD->NumEle) + i*MD->NumRiv + j] = Y[(3*MD->NumEle+2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j] / (Y[3*MD->NumEle+j]*MD->Riv[j%MD->NumRiv].coeff*MD->Riv[j%MD->NumRiv].Length);
		     //DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = Y[(3*MD->NumEle+2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]; //?? TODO
		    DummyAgeC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = MD->RivAgeConc[j][i];
	     else
		     DummyAgeC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = 0.0;
		//printf(" %lf\t",DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j]);
		if(j==0*MD->NumRiv) printf(" %lf\t",DummyAgeC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j]);
        }
    }
	//??
	/*
    if(t<=10){
    	DummyC[0 + 2*MD->NumEle+667-1] = 1.0; //??
	DummyC[667-1] = 2.0; //??
	DummyC[MD->NumSolute*3*MD->NumEle+0 + 0*MD->NumRiv] = 3.0; //??
    }
    else{
	DummyC[0 + 2*MD->NumEle+667-1] = 0.0; //??
        DummyC[667-1] = 0.0; //??
        DummyC[MD->NumSolute*3*MD->NumEle+0 + 0*MD->NumRiv] = 0.0; //??
    }
	//??*/

}

/* AGE: Initialize AgeC terms - END */






int f_sol(FILE **omf, realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
	{




	// TRANSPORT VARIABLES - START /
	realtype *DummyC, *DummyDC, *DummyAgeC, *DummyAgeDC;
	realtype CPrep, CInfil, AgeCPrep, AgeCInfil;
	realtype **uSub, **uSurf, **qRiv, **uRiv, *uRivCount, uShear, HydraulicRadius, *uPrep, *uInfilUsat, *uInfilSat, *uETSurf, *uETUsat, *uETSat;
	realtype Cs, As, Ab, Ab_down;
    	realtype Ux, Uy, U, uCount, AlphaL, AlphaT, Diffusivity;
    	realtype **Disp, *DispRiv;
    	realtype Dxx, Dxy, Dyx, Dyy, Cos, Sin;
    	realtype Grad_C, Dist;
    	realtype AgeCs, Grad_AgeC;
	realtype tDC;
	realtype **Kd, **R;
	// TRANSPORT VARIABLES - END /
// Yi
        realtype UsatDepth;
        realtype *DummyS,*DummyDS,*DummySmas,*DummyDSmas;
        int cid;
        realtype tttt, tttt1,tttt2,tttt3;
        realtype somsobT,somsolT,xmas;
// Yi

  	int i, j, k, inabr;
  	realtype Delta, Gamma;
  	  	realtype Rn, G, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a,r_s,alpha_r,f_r,eta_s,beta_s,gamma_s,Rmax, Lambda, P_c,qv, qv_sat, ETp;
	realtype ThetaRef, ThetaW;
  	realtype Avg_Y_Surf, Dif_Y_Surf,Grad_Y_Surf, Avg_Sf,Distance;
  	realtype Cwr,TotalY_Riv, TotalY_Riv_down,CrossA,CrossAdown,AvgCrossA,Perem, Perem_down,Avg_Rough,Avg_Perem,Avg_Y_Riv,Dif_Y_Riv,Grad_Y_Riv,Wid,Wid_down,Avg_Wid;
  	realtype Avg_Y_Sub, Dif_Y_Sub,Avg_Ksat, Grad_Y_Sub,nabrAqDepth,AquiferDepth, Deficit,elemSatn,satKfunc,effK,effKnabr,TotalY_Ele,TotalY_Ele_down;
  	realtype *Y, *DY;
  	Model_Data MD;
  	Y = NV_DATA_S(CV_Y);
  	//DY = NV_DATA_S(CV_Ydot);
  	MD = (Model_Data) DS;
	//TODO Fix 10000
	DY = (realtype *)malloc(10000*sizeof(realtype));
	//printf("\n$");
	for(i=0; i<MD->NumRiv; i++)
		//printf("%lf\t", Y[6*MD->NumEle+2*MD->NumRiv+i]);
		//printf("%lf\t",NV_Ith_S(CV_Y,6*MD->NumEle+2*MD->NumRiv+i));

	//?#printf("\nHello from f_sol.c\n"); getchar();
	/* Initialization of temporary state variables */
	for(i=0; i<3*MD->NumEle+2*MD->NumRiv; i++)
  		{
		MD->DummyY[i]=(Y[i]>=0)?Y[i]:0;
		DY[i]=0;
  		if(i<MD->NumRiv)
			{
			MD->FluxRiv[i][0]=0;
			MD->FluxRiv[i][10]=0;
			}
		if((MD->SurfMode==2)&&(i<MD->NumEle))
			{
			for(j=0;j<3;j++)
				{
				MD->Ele[i].surfH[j]=(MD->Ele[i].nabr[j]>0)?((MD->Ele[i].BC[j]>-4)?(MD->Ele[MD->Ele[i].nabr[j]-1].zmax+MD->DummyY[MD->Ele[i].nabr[j]-1]):((MD->DummyY[-(MD->Ele[i].BC[j]/4)-1+3*MD->NumEle]>MD->Riv[-(MD->Ele[i].BC[j]/4)-1].depth)?MD->Riv[-(MD->Ele[i].BC[j]/4)-1].zmin+MD->DummyY[-(MD->Ele[i].BC[j]/4)-1+3*MD->NumEle]:MD->Riv[-(MD->Ele[i].BC[j]/4)-1].zmax)):((MD->Ele[i].BC[j]!=1)?(MD->Ele[i].zmax+MD->DummyY[i]):Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t));
				}
                	MD->Ele[i].dhBYdx=-1*(MD->Ele[i].surfY[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfX[2]*(MD->Ele[i].surfY[1]-MD->Ele[i].surfY[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfY[0]-MD->Ele[i].surfY[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfY[2]-MD->Ele[i].surfY[1]));
                	MD->Ele[i].dhBYdy=-1*(MD->Ele[i].surfX[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfY[2]*(MD->Ele[i].surfX[1]-MD->Ele[i].surfX[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfX[0]-MD->Ele[i].surfX[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfX[2]-MD->Ele[i].surfX[1]));
			}
  		}	






/* TRANSPORT - INITIALIZE TRANSPORT VARIABLES - START */
    DummyC = (realtype *)malloc((3*MD->NumEle+2*MD->NumRiv)*MD->NumSolute*sizeof(realtype));
    DummyDC= (realtype *)malloc((3*MD->NumEle+2*MD->NumRiv)*MD->NumSolute*sizeof(realtype));
    DummyAgeC = (realtype *)malloc((3*MD->NumEle+2*MD->NumRiv)*MD->NumSolute*sizeof(realtype));
    DummyAgeDC= (realtype *)malloc((3*MD->NumEle+2*MD->NumRiv)*MD->NumSolute*sizeof(realtype));
//Yi
    DummyS = (realtype *)malloc((2*MD->NumEle)*MD->NumSolute*sizeof(realtype));
    DummyDS = (realtype *)malloc((2*MD->NumEle)*MD->NumSolute*sizeof(realtype));
    DummySmas = (realtype *)malloc((2*MD->NumEle)*MD->NumSolute*sizeof(realtype));
    DummyDSmas = (realtype *)malloc((2*MD->NumEle)*MD->NumSolute*sizeof(realtype));
//Yi
    uSub  = (realtype **)malloc(MD->NumEle*sizeof(realtype *));
    uSurf = (realtype **)malloc(MD->NumEle*sizeof(realtype *));
    qRiv  = (realtype **)malloc(MD->NumRiv*sizeof(realtype *));
    uRiv  = (realtype **)malloc(MD->NumRiv*sizeof(realtype *));
    uRivCount = (realtype *)malloc(MD->NumRiv*sizeof(realtype));

    uPrep   = (realtype *)malloc(MD->NumEle*sizeof(realtype));
    uInfilUsat  = (realtype *)malloc(MD->NumEle*sizeof(realtype));
    uInfilSat   = (realtype *)malloc(MD->NumEle*sizeof(realtype));
    uETSurf = (realtype *)malloc(MD->NumEle*sizeof(realtype));
    uETUsat  = (realtype *)malloc(MD->NumEle*sizeof(realtype));
    uETSat  = (realtype *)malloc(MD->NumEle*sizeof(realtype));

    for(i=0; i<MD->NumEle; i++){
	uSub[i] = (realtype *)malloc(3*sizeof(realtype));
    	uSurf[i]= (realtype *)malloc(3*sizeof(realtype));
    }
    for(i=0; i<MD->NumRiv; i++){
	//qRiv[i] = (realtype *)malloc(10*sizeof(realtype));
	qRiv[i] = (realtype *)malloc(11*sizeof(realtype)); //Yi
	uRiv[i] = (realtype *)malloc(2*sizeof(realtype));
	uRivCount[i] = 0;
    }

    Disp = (realtype **)malloc(MD->NumEle*sizeof(realtype *));
    for(i=0; i<MD->NumEle; i++)
	Disp[i] = (realtype *)malloc(3*sizeof(realtype));
    DispRiv = (realtype *)malloc(MD->NumRiv*sizeof(realtype));


    Kd = (realtype **)malloc(MD->NumSolute*sizeof(realtype *));
    R  = (realtype **)malloc(MD->NumSolute*sizeof(realtype *));
    for(i=0; i<MD->NumSolute; i++){
	Kd[i] = (realtype *)malloc(MD->NumEle*sizeof(realtype));
	R[i]  = (realtype *)malloc(MD->NumEle*sizeof(realtype));
    }
    

	//??InitFlag = checkReInit(MD, CV_Y, Y); if (InitFlag == 1) return (1);
	realtype C;
	C = Interpolation(&MD->TSD_PChem[0*MD->NumPrep+MD->Ele[0].prep-1], t); // TODO Fix ME 
	//printf("%lf", C); getchar();
    initDummyC(DummyC, MD, Y, t, C); AGE initDummyAgeC( DummyAgeC,  MD, Y, t);
    initDummyDC(DummyDC, MD);     AGE initDummyAgeDC(DummyAgeDC, MD);
//Yi

    for (i = 0; i< MD->NumEle; i++){
            AquiferDepth=(MD->Ele[i].zmax-MD->Ele[i].zmin);
            DummyS[i + 0*MD->NumEle] = MD->EleSbar[i + 0*MD->NumEle];
            DummyS[i + 1*MD->NumEle] = MD->EleSbar[i + 1*MD->NumEle];
            DummySmas[i + 0*MD->NumEle] = MD->EleS[i + 0*MD->NumEle];
            DummySmas[i + 1*MD->NumEle] = MD->EleS[i + 1*MD->NumEle];

            DummyDS[i + 0*MD->NumEle] = 0.0;
            DummyDS[i + 1*MD->NumEle] = 0.0;
	    DummyDSmas[i + 0*MD->NumEle] = 0.0;
            DummyDSmas[i + 1*MD->NumEle] = 0.0;
    } 
   cid = 53;
//Yi
/* TRANSPORT - INITIALIZE TRANSPORT VARIABLES - END */









/* Lateral Flux Calculation between Triangular elements Follows  */
	for(i=0; i<MD->NumEle; i++)
  		{
		AquiferDepth=(MD->Ele[i].zmax-MD->Ele[i].zmin);
		if(AquiferDepth<MD->Ele[i].macD) MD->Ele[i].macD=AquiferDepth;
    		for(j=0; j<3; j++)
    			{
      			if(MD->Ele[i].nabr[j] > 0)
      				{
				/***************************************************************************/
				/* Subsurface Lateral Flux Calculation between Triangular elements Follows */
				/***************************************************************************/
        			Dif_Y_Sub = (MD->DummyY[i+2*MD->NumEle] + MD->Ele[i].zmin) - (MD->DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle] + MD->Ele[MD->Ele[i].nabr[j]-1].zmin);
//				Avg_Y_Sub=avgY(MD->Ele[i].zmin,MD->Ele[MD->Ele[i].nabr[j]-1].zmin,MD->DummyY[i+2*MD->NumEle],MD->DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle]);
				Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[i+2*MD->NumEle],MD->DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle]);
        			Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[MD->Ele[i].nabr[j] - 1].x), 2) + pow((MD->Ele[i].y - MD->Ele[MD->Ele[i].nabr[j] - 1].y), 2));
        			Grad_Y_Sub = Dif_Y_Sub/Distance;
        			/* take care of macropore effect */
				effK=effKH(MD->Ele[i].Macropore,MD->DummyY[i+2*MD->NumEle],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH);
				inabr=MD->Ele[i].nabr[j]-1;
				nabrAqDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
				effKnabr=effKH(MD->Ele[inabr].Macropore,MD->DummyY[inabr+2*MD->NumEle],nabrAqDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
				/* It should be weighted average. However, there is an ambiguity about distance used */
         			Avg_Ksat=0.5*(effK+effKnabr); 
        			/* groundwater flow modeled by Darcy's law */
        			MD->FluxSub[i][j] = Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub*MD->Ele[i].edge[j];        
				/***************************************************************************/
				/* Surface Lateral Flux Calculation between Triangular elements Follows    */
				/***************************************************************************/
        			Dif_Y_Surf = (MD->SurfMode==1)?(MD->Ele[i].zmax-MD->Ele[MD->Ele[i].nabr[j] - 1].zmax):(MD->DummyY[i] + MD->Ele[i].zmax) - (MD->DummyY[MD->Ele[i].nabr[j] - 1] + MD->Ele[MD->Ele[i].nabr[j] - 1].zmax);
//				Avg_Y_Surf=avgY(MD->Ele[i].zmax,MD->Ele[MD->Ele[i].nabr[j] - 1].zmax,MD->DummyY[i],MD->DummyY[MD->Ele[i].nabr[j]-1]);
				Avg_Y_Surf=avgY(Dif_Y_Surf,MD->DummyY[i],MD->DummyY[MD->Ele[i].nabr[j]-1]);
        			Grad_Y_Surf = Dif_Y_Surf/Distance;
				Avg_Sf=0.5*(sqrt(pow(MD->Ele[i].dhBYdx,2)+pow(MD->Ele[i].dhBYdy,2)) + sqrt(pow(MD->Ele[MD->Ele[i].nabr[j] - 1].dhBYdx, 2) + pow(MD->Ele[MD->Ele[i].nabr[j] - 1].dhBYdy, 2)));//Yi weight
				Avg_Sf=(MD->SurfMode==1)?(Grad_Y_Surf>0?Grad_Y_Surf:EPS/pow(10.0,6)):(Avg_Sf>EPS/pow(10.0,6))?Avg_Sf:EPS/pow(10.0,6);
				/* Weighting needed */
        			Avg_Rough = 0.5*(MD->Ele[i].Rough + MD->Ele[MD->Ele[i].nabr[j] - 1].Rough);
        			CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];        
				OverlandFlow(MD->FluxSurf,i,j, Avg_Y_Surf,Grad_Y_Surf,Avg_Sf,CrossA,Avg_Rough);
				//??MD->FluxSurf[i][j]=0.0;//??bhatt
      				}
			/************************************************/
			/* Boundary condition Flux Calculations Follows */
			/************************************************/
      			else
      				{
        			/*  No flow (natural) boundary condition is default */
        			if(MD->Ele[i].BC[j] == 0)
        				{
          				MD->FluxSurf[i][j] = 0;
          				MD->FluxSub[i][j] = 0;
        				}
        			else if(MD->Ele[i].BC[j] == 1)	/* Note: ideally different boundary conditions need to be incorporated	for surf and subsurf respectively */
					/* Note: the formulation assumes only dirichlet TS right now */
        				{
          				MD->FluxSurf[i][j] = 0;	/* Note the assumption here is no flow for surface*/ 
            				Dif_Y_Sub = (MD->DummyY[i+2*MD->NumEle] + MD->Ele[i].zmin) - Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t);            
            				Avg_Y_Sub = avgY(Dif_Y_Sub,MD->DummyY[i+2*MD->NumEle],(Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t) - MD->Ele[i].zmin));
//            				Avg_Y_Sub = (MD->DummyY[i+2*MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t) - MD->Ele[i].zmin))/2;
	    				/* Minimum Distance from circumcenter to the edge of the triangle on which BDD. condition is defined*/
            				Distance = sqrt(pow(MD->Ele[i].edge[0]*MD->Ele[i].edge[1]*MD->Ele[i].edge[2]/(4*MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j]/2, 2));
					effK=effKH(MD->Ele[i].Macropore,MD->DummyY[i+2*MD->NumEle],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH);
            				Avg_Ksat = effK;
            				Grad_Y_Sub = Dif_Y_Sub/Distance;
            				MD->FluxSub[i][j] = Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub*MD->Ele[i].edge[j];
          				}
          			else       /* Neumann BC (Note: MD->Ele[i].BC[j] value have to be = 2+(index of neumann boundary TS)*/
          				{
            				MD->FluxSurf[i][j] = Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t);
            				MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC[j])-1], t);
          				}
      				}
    			} 
		/**************************************************************************************************/
		/* Evaporation Module: [2] is ET from OVLF/SUBF, [1] is Transpiration, [0] is ET loss from canopy */
		/**************************************************************************************************/
		/* Physical Unit Dependent. Change this */
    		Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn-1], t);
	//    	G = Interpolation(&MD->TSD_G[MD->Ele[i].G-1], t);
			G = 0.1*Rn;
    		T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp-1], t);
    		Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel-1], t);
    		RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity-1], t);
    				VP = 611.2*exp(17.67*T/(T+243.5))*RH;
    		P = 101.325*pow(10,3)*pow((293-0.0065*MD->Ele[i].zmax)/293,5.26);
		qv = 0.622*VP/P;
		qv_sat = 0.622*(VP/RH)/P;
    		//P = 101.325*pow(10,3)*pow((293-0.0065*MD->Ele[i].zmax)/293,5.26);
    		//Delta = 2503*pow(10,3)*exp(17.27*T/(T+237.3))/(pow(237.3 + T, 2));
    		//Gamma = P*1.0035*0.92/(0.622*2441);
    		LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC-1], t);
/*    		zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
    		cnpy_h = zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25))));
    		if(LAI<2.85)	
    			{
    			rl= 0.0002 + 0.3*cnpy_h*pow(0.07*LAI,0.5);
    			}
    		else
    			{
    			rl= 0.3*cnpy_h*(1-(zero_dh/cnpy_h));
    			} */
  		rl=Interpolation(&MD->TSD_RL[MD->Ele[i].LC-1], t);
    		r_a = 12*4.72*log(MD->Ele[i].windH/rl)/(0.54*Vel/UNIT_C/60+1)/UNIT_C/60;

		Gamma = 4*0.7*SIGMA*UNIT_C*R_dry/C_air*pow(T+273.15,4)/(P/r_a)+1;
		Delta = Lv*Lv*0.622/R_v/C_air/pow(T+273.15,2)*qv_sat;
		ETp = (Rn*Delta+Gamma*(1.2*Lv*(qv_sat-qv)/r_a))/(1000.0*Lv*(Delta+Gamma));   
    		//MD->EleET[i][2] = MD->pcCal.Et2*(1-MD->Ele[i].VegFrac)*(Rn*(1-MD->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a))/(1000.0*2441000.0*(Delta+Gamma));
   		if((MD->Ele[i].zmax-MD->Ele[i].zmin)-MD->DummyY[i+2*MD->NumEle]<MD->Ele[i].RzD)
   			{
   			elemSatn=1.0;	
   			}
   		else
   			{
			elemSatn = ((MD->DummyY[i+MD->NumEle]/(AquiferDepth-MD->DummyY[i+2*MD->NumEle]))>1)?1:((MD->DummyY[i+MD->NumEle]/(AquiferDepth-MD->DummyY[i+2*MD->NumEle]))<0)?0:0.5*(1-cos(3.14*(MD->DummyY[i+MD->NumEle]/(AquiferDepth-MD->DummyY[i+2*MD->NumEle])))); 
   			}
   			ThetaRef = 0.7*MD->Soil[(MD->Ele[i].soil-1)].ThetaS;
   				ThetaW = 1.05*MD->Soil[(MD->Ele[i].soil-1)].ThetaR;
                beta_s= (elemSatn*MD->Ele[i].Porosity+MD->Soil[(MD->Ele[i].soil-1)].ThetaR-ThetaW)/(ThetaRef-ThetaW);
                beta_s = (beta_s<0.0001)?0.0001:(beta_s>1?1:beta_s);
                MD->EleET[i][2] = MD->pcCal.Et2*(1-MD->Ele[i].VegFrac)*beta_s*ETp;
		MD->EleET[i][2] = MD->EleET[i][2]<0?0:MD->EleET[i][2];
    		if(LAI>0.0)
    			{
    	    		Rmax = 5000.0/(60*UNIT_C);		/* Unit day_per_m */
	    			f_r= 1.1*1.5*Rn/(MD->Ele[i].Rs_ref*LAI);
			f_r = f_r<0?0:f_r;
	    		alpha_r= (1+f_r)/(f_r+(MD->Ele[i].Rmin/Rmax));
			alpha_r = alpha_r>10000?10000:alpha_r;
	    		eta_s= 1- 0.0016*(pow((24.85-T),2));
			eta_s=eta_s<0.0001?0.0001:eta_s;
			gamma_s=1/(1+0.00025*(VP/RH-VP));
			gamma_s = (gamma_s<0.01)?0.01:gamma_s;
	    		r_s=((MD->Ele[i].Rmin*alpha_r/(beta_s*LAI*eta_s*gamma_s))> Rmax)?Rmax:(MD->Ele[i].Rmin*alpha_r/(beta_s*LAI*eta_s*gamma_s));
			P_c = (1+Delta/Gamma)/(1+r_s/r_a+Delta/Gamma);
	    		MD->EleET[i][1] = MD->pcCal.Et1*MD->Ele[i].VegFrac*P_c*(1-pow(((MD->EleIS[i]+MD->EleSnowCanopy[i]<0)?0:(MD->EleIS[i]+MD->EleSnowCanopy[i]))/(MD->EleISmax[i]+MD->EleISsnowmax[i]),1.0/2.0))*ETp;
			MD->EleET[i][1] = MD->EleET[i][1]<0?0:MD->EleET[i][1];
    			AquiferDepth=MD->Ele[i].zmax-MD->Ele[i].zmin;//??BHATT
                        MD->EleET[i][1] = ((MD->DummyY[i+2*MD->NumEle]<(AquiferDepth-MD->Ele[i].RzD))&&MD->DummyY[i+MD->NumEle]<=0)?0:MD->EleET[i][1];//??BHATT
			}
    		else
    			{
    			MD->EleET[i][1] =0.0;
    			}
		/* Note: Assumption is OVL flow depth less than EPS/100 is immobile water */
		if(MD->DummyY[i+2*MD->NumEle]>AquiferDepth-MD->Ele[i].infD)
      			{
			/* Assumption: infD<macD */
			Grad_Y_Sub=(MD->DummyY[i]+MD->Ele[i].zmax-(MD->DummyY[i+2*MD->NumEle]+MD->Ele[i].zmin))/MD->Ele[i].infD;
			Grad_Y_Sub=((MD->DummyY[i]<EPS/100)&&(Grad_Y_Sub>0))?0:Grad_Y_Sub;
			elemSatn=1.0;
			satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,MD->Ele[i].Beta/(MD->Ele[i].Beta-1)),(MD->Ele[i].Beta-1)/MD->Ele[i].Beta),2);
			effK=(1)?effKV(satKfunc,Grad_Y_Sub,MD->Ele[i].macKsatV,MD->Ele[i].infKsatV,MD->Ele[i].hAreaF):MD->Ele[i].infKsatV;
      			MD->EleViR[i] = effK*Grad_Y_Sub;
			MD->Recharge[i] = MD->EleViR[i];
       			DY[i+MD->NumEle] = DY[i+MD->NumEle]+MD->EleViR[i]-MD->Recharge[i];
  			DY[i+2*MD->NumEle]=DY[i+2*MD->NumEle]+MD->Recharge[i]-((MD->DummyY[i]<EPS/100)?MD->EleET[i][2]:0);
//??			if(DY[i+2*MD->NumEle]<0 && MD->DummyY[i+2*MD->NumEle]<.002) printf("1 %d %lf %lf\n", MD->Ele[i].soil, MD->Recharge[i], MD->EleET[i][2]); //??BHATT
			}
      		else
      			{
			Deficit=AquiferDepth-MD->DummyY[i+2*MD->NumEle];
			elemSatn=((MD->DummyY[i+MD->NumEle]/Deficit)>1)?1:((MD->DummyY[i+MD->NumEle]<=0)?EPS/1000.0:MD->DummyY[i+MD->NumEle]/Deficit);
			/* Note: for psi calculation using van genuchten relation, cutting the psi-sat tail at small saturation can be performed for computational advantage. If you dont' want to perform this, comment the statement that follows */
			elemSatn=(elemSatn<multF*EPS)?multF*EPS:elemSatn;
			Avg_Y_Sub=(-(pow(pow(1/elemSatn,MD->Ele[i].Beta/(MD->Ele[i].Beta-1))-1,1/MD->Ele[i].Beta)/MD->Ele[i].Alpha)<MINpsi)?MINpsi:-(pow(pow(1/elemSatn,MD->Ele[i].Beta/(MD->Ele[i].Beta-1))-1,1/MD->Ele[i].Beta)/MD->Ele[i].Alpha);
			TotalY_Ele=Avg_Y_Sub+MD->Ele[i].zmin+AquiferDepth-MD->Ele[i].infD;
			Grad_Y_Sub=(MD->DummyY[i]+MD->Ele[i].zmax-TotalY_Ele)/MD->Ele[i].infD;
			Grad_Y_Sub=((MD->DummyY[i]<EPS/100)&&(Grad_Y_Sub>0))?0:Grad_Y_Sub;
			satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,MD->Ele[i].Beta/(MD->Ele[i].Beta-1)),(MD->Ele[i].Beta-1)/MD->Ele[i].Beta),2);
			satKfunc=(satKfunc<0.13)?0.13:satKfunc;//Xuan??
			effK=(1)?effKV(satKfunc,Grad_Y_Sub,MD->Ele[i].macKsatV,MD->Ele[i].infKsatV,MD->Ele[i].hAreaF):MD->Ele[i].infKsatV;
     	//		MD->EleViR[i] = 0.5*(effK+MD->Ele[i].infKsatV)*Grad_Y_Sub;//BHATT??
		MD->EleViR[i] = (effK)*Grad_Y_Sub;
			/* Harmonic mean formulation. Note that if unsaturated zone has low saturation, satKfunc becomes very small. Use arithmetic mean instead*/
//	                MD->Recharge[i] = (elemSatn==0.0)?0:(Deficit<=0)?0:(MD->Ele[i].KsatV*satKfunc*(MD->Ele[i].Alpha*Deficit-2*pow(-1+pow(elemSatn,MD->Ele[i].Beta/(-MD->Ele[i].Beta+1)),1/MD->Ele[i].Beta))/(MD->Ele[i].Alpha*((Deficit+MD->DummyY[i+2*MD->NumEle]*satKfunc))));
			/* Arithmetic Mean Formulation */
			effK=(MD->Ele[i].Macropore==1)?((MD->DummyY[i+2*MD->NumEle]>AquiferDepth-MD->Ele[i].macD)?effK:MD->Ele[i].KsatV*satKfunc):MD->Ele[i].KsatV*satKfunc;
	                MD->Recharge[i] = (elemSatn==0.0)?0:(Deficit<=0)?0:(MD->Ele[i].KsatV*MD->DummyY[i+2*MD->NumEle]+effK*Deficit)*(MD->Ele[i].Alpha*Deficit-2*pow(-1+pow(elemSatn,MD->Ele[i].Beta/(-MD->Ele[i].Beta+1)),1/MD->Ele[i].Beta))/(MD->Ele[i].Alpha*pow(Deficit+MD->DummyY[i+2*MD->NumEle],2));
			MD->Recharge[i]=(MD->Recharge[i]>0&&MD->DummyY[i+MD->NumEle]<=0)?0:MD->Recharge[i]; //?? BHATT
			MD->Recharge[i]=(MD->Recharge[i]<0&&MD->DummyY[i+2*MD->NumEle]<=0)?0:MD->Recharge[i];//?? BHATT
			MD->EleET[i][2]=(MD->DummyY[i]<EPS/100)?elemSatn*MD->EleET[i][2]:MD->EleET[i][2];
			MD->EleET[i][2]=(MD->DummyY[i]<EPS/100)?(elemSatn<=multF*EPS?0:MD->EleET[i][2]):MD->EleET[i][2];//??BHATT
       			DY[i+MD->NumEle] = DY[i+MD->NumEle]+MD->EleViR[i]-MD->Recharge[i]-((MD->DummyY[i]<EPS/100)?MD->EleET[i][2]:0);	
  			DY[i+2*MD->NumEle]=DY[i+2*MD->NumEle]+MD->Recharge[i];
			}
       		DY[i] = DY[i]+MD->EleNetPrep[i] - MD->EleViR[i]-((MD->DummyY[i]<EPS/100)?0:MD->EleET[i][2]);
		if(MD->DummyY[i+2*MD->NumEle]>AquiferDepth-MD->Ele[i].RzD)
			{
			DY[i+2*MD->NumEle]=DY[i+2*MD->NumEle]-MD->EleET[i][1];
			}
		else
			{
			DY[i+MD->NumEle] = DY[i+MD->NumEle]-MD->EleET[i][1];
			}
//??		if(DY[i+MD->NumEle]<0 && MD->DummyY[i+MD->NumEle]<.002) printf("2 %d %lf %lf %lf %lf\n", MD->Ele[i].soil, MD->EleViR[i],MD->Recharge[i], MD->EleET[i][2], MD->EleET[i][1]); //??BHATT
//??            if(DY[i+2*MD->NumEle]<0 && MD->DummyY[i+2*MD->NumEle]<.002) {printf("3 %d %lf %lf\n", MD->Ele[i].soil, MD->Recharge[i], MD->EleET[i][1]); getchar(); }//?? BHATT
  		}
	/* Lateral Flux Calculation between River-River and River-Triangular elements Follows */   
  	for(i=0; i<MD->NumRiv; i++)
  		{
    		TotalY_Riv = MD->DummyY[i + 3*MD->NumEle] + MD->Riv[i].zmin;
    		Perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,2);
    		if(MD->Riv[i].down > 0)
    			{
			/****************************************************************/
			/* Lateral Flux Calculation between River-River element Follows */    
			/****************************************************************/
      			TotalY_Riv_down = MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle] + MD->Riv[MD->Riv[i].down - 1].zmin;
      			Perem_down = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd,MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle],MD->Riv[MD->Riv[i].down - 1].coeff,2);
      			Avg_Perem = (Perem + Perem_down)/2.0;	
      		//	Avg_Rough = (MD->Riv_Mat[MD->Riv[i].material - 1].Rough + MD->Riv_Mat[MD->Riv[MD->Riv[i].down - 1].material-1].Rough)/2.0;
                        Avg_Rough = (MD->Riv[i].Rough + MD->Riv[MD->Riv[i].down - 1].Rough) / 2.0; //Yi
      			Distance = 0.5*(MD->Riv[i].Length+MD->Riv[MD->Riv[i].down - 1].Length);
			Dif_Y_Riv=(MD->RivMode==1)?(MD->Riv[i].zmin-MD->Riv[MD->Riv[i].down - 1].zmin):(TotalY_Riv - TotalY_Riv_down);
      			Grad_Y_Riv = Dif_Y_Riv/Distance; 
      			Avg_Sf = (Grad_Y_Riv>0)?Grad_Y_Riv:EPS;   
			CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,1);      
			CrossAdown = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd,MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle],MD->Riv[MD->Riv[i].down - 1].coeff,1);     
			AvgCrossA=0.5*(CrossA+CrossAdown); 
			Avg_Y_Riv=(Avg_Perem==0)?0:(AvgCrossA/Avg_Perem);
			OverlandFlow(MD->FluxRiv,i,1, Avg_Y_Riv,Grad_Y_Riv,Avg_Sf,CrossA,Avg_Rough);
      			/* accumulate to get in-flow for down segments: [0] for inflow, [1] for outflow */
      			MD->FluxRiv[MD->Riv[i].down - 1][0] = MD->FluxRiv[MD->Riv[i].down - 1][0] - MD->FluxRiv[i][1];        
			/************************************************************************/
			/* Lateral Flux Calculation between Element Beneath River (EBR) and EBR */    
			/************************************************************************/
      			TotalY_Ele = MD->DummyY[i + 3*MD->NumEle+MD->NumRiv] + MD->Ele[i+MD->NumEle].zmin;
      			TotalY_Ele_down = MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle+MD->NumRiv] + MD->Ele[MD->Riv[i].down - 1+MD->NumEle].zmin;
    			Wid = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3);
      			Wid_down = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd,MD->Riv[MD->Riv[i].down - 1].depth,MD->Riv[MD->Riv[i].down - 1].coeff,3);
      			Avg_Wid = (Wid + Wid_down)/2.0;	
      			Distance = 0.5*(MD->Riv[i].Length+MD->Riv[MD->Riv[i].down - 1].Length);
			Dif_Y_Sub=TotalY_Ele - TotalY_Ele_down;
//     			Avg_Y_Sub=avgY(MD->Ele[i+MD->NumEle].zmin,MD->Ele[MD->Riv[i].down - 1+MD->NumEle].zmin,MD->DummyY[i + 3*MD->NumEle+MD->NumRiv],MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle+MD->NumRiv]);
      			Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[i + 3*MD->NumEle+MD->NumRiv],MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle+MD->NumRiv]);
      			Grad_Y_Sub = Dif_Y_Sub/Distance; 
                       	/* take care of macropore effect */
			AquiferDepth=MD->Ele[i+MD->NumEle].zmax-MD->Ele[i+MD->NumEle].zmin;
//                     	effK=MD->Ele[i+MD->NumEle].KsatH;
                     	effK=0.5*(effKH(MD->Ele[MD->Riv[i].LeftEle-1].Macropore,MD->DummyY[MD->Riv[i].LeftEle-1+2*MD->NumEle],MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin,MD->Ele[MD->Riv[i].LeftEle-1].macD,MD->Ele[MD->Riv[i].LeftEle-1].macKsatH,MD->Ele[MD->Riv[i].LeftEle-1].vAreaF,MD->Ele[MD->Riv[i].LeftEle-1].KsatH)+effKH(MD->Ele[MD->Riv[i].RightEle-1].Macropore,MD->DummyY[MD->Riv[i].RightEle-1+2*MD->NumEle],MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin,MD->Ele[MD->Riv[i].RightEle-1].macD,MD->Ele[MD->Riv[i].RightEle-1].macKsatH,MD->Ele[MD->Riv[i].RightEle-1].vAreaF,MD->Ele[MD->Riv[i].RightEle-1].KsatH));
                     	inabr=MD->Riv[i].down - 1;
                        nabrAqDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
//                     	effKnabr=MD->Ele[inabr+MD->NumEle].KsatH;
                     	effKnabr=0.5*(effKH(MD->Ele[MD->Riv[inabr].LeftEle-1].Macropore,MD->DummyY[MD->Riv[inabr].LeftEle-1+2*MD->NumEle],MD->Ele[MD->Riv[inabr].LeftEle-1].zmax-MD->Ele[MD->Riv[inabr].LeftEle-1].zmin,MD->Ele[MD->Riv[inabr].LeftEle-1].macD,MD->Ele[MD->Riv[inabr].LeftEle-1].macKsatH,MD->Ele[MD->Riv[inabr].LeftEle-1].vAreaF,MD->Ele[MD->Riv[inabr].LeftEle-1].KsatH)+effKH(MD->Ele[MD->Riv[inabr].RightEle-1].Macropore,MD->DummyY[MD->Riv[inabr].RightEle-1+2*MD->NumEle],MD->Ele[MD->Riv[inabr].RightEle-1].zmax-MD->Ele[MD->Riv[inabr].RightEle-1].zmin,MD->Ele[MD->Riv[inabr].RightEle-1].macD,MD->Ele[MD->Riv[inabr].RightEle-1].macKsatH,MD->Ele[MD->Riv[inabr].RightEle-1].vAreaF,MD->Ele[MD->Riv[inabr].RightEle-1].KsatH));
                    	Avg_Ksat=0.5*(effK+effKnabr);
                     	/* groundwater flow modeled by Darcy's law */
                    	MD->FluxRiv[i][9] = Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub*Avg_Wid;
      			/* accumulate to get in-flow for down segments: [10] for inflow, [9] for outflow */
      			MD->FluxRiv[MD->Riv[i].down - 1][10] = MD->FluxRiv[MD->Riv[i].down - 1][10] - MD->FluxRiv[i][9];        
    			} 
    		else
    			{
      			switch(MD->Riv[i].down)
      				{
        			case -1:  
          			/* Dirichlet boundary condition */        
          			TotalY_Riv_down = Interpolation(&MD->TSD_Riv[(MD->Riv[i].BC)-1], t) + (MD->Node[MD->Riv[i].ToNode-1].zmax -MD->Riv[i].depth);
          			Distance = sqrt(pow(MD->Riv[i].x - MD->Node[MD->Riv[i].ToNode-1].x, 2) + pow(MD->Riv[i].y - MD->Node[MD->Riv[i].ToNode-1].y, 2));
          			Grad_Y_Riv = (TotalY_Riv - TotalY_Riv_down)/Distance; 
				/* Note: do i need to change else part here for diff wave */
      				Avg_Sf = (MD->RivMode==1)?Grad_Y_Riv:Grad_Y_Riv;;   
          			Avg_Rough = MD->Riv_Mat[MD->Riv[i].material-1].Rough;
                                Avg_Rough = MD->Riv[i].Rough; //Yi
          			Avg_Y_Riv = avgY(Grad_Y_Riv,MD->DummyY[i + 3*MD->NumEle],Interpolation(&MD->TSD_Riv[(MD->Riv[i].BC)-1], t));
          			Avg_Perem = Perem;
          			CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,1);
                        	Avg_Y_Riv=(Perem==0)?0:(CrossA/Avg_Perem);
                        	OverlandFlow(MD->FluxRiv,i,1, Avg_Y_Riv,Grad_Y_Riv,Avg_Sf,CrossA,Avg_Rough);			
          			break;
        			case -2:     
          			/* Neumann boundary condition */
          			MD->FluxRiv[i][1] = Interpolation(&MD->TSD_Riv[MD->Riv[i].BC-1], t);
          			break;  
        			case -3:     
          			/* zero-depth-gradient boundary conditions */
          			Distance = sqrt(pow(MD->Riv[i].x - MD->Node[MD->Riv[i].ToNode-1].x, 2) + pow(MD->Riv[i].y - MD->Node[MD->Riv[i].ToNode-1].y, 2));
          			Grad_Y_Riv = (MD->Riv[i].zmin - (MD->Node[MD->Riv[i].ToNode-1].zmax -MD->Riv[i].depth))/Distance;
          			Avg_Rough = MD->Riv_Mat[MD->Riv[i].material-1].Rough;
                                Avg_Rough = MD->Riv[i].Rough; //Yi
          			Avg_Y_Riv = MD->DummyY[i + 3*MD->NumEle];
          			Avg_Perem = Perem;
          			CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,1);
          			MD->FluxRiv[i][1] = sqrt(Grad_Y_Riv)*CrossA*((Avg_Perem>0)?pow(CrossA/Avg_Perem,2.0/3.0):0)/Avg_Rough; 
          			break;
        			case -4:          
          			/* Critical Depth boundary conditions */
          			CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,1);
          			MD->FluxRiv[i][1] = CrossA*sqrt(GRAV*UNIT_C*UNIT_C*MD->DummyY[i + 3*MD->NumEle]);  /* Note the dependence on physical units */
         	 		break;
        			default:
          			printf("Fatal Error: River Routing Boundary Condition Type Is Wrong!");
          			exit(1); 
      				}   
			/* Note: bdd condition for subsurface element can be changed. Assumption: No flow condition */ 
			MD->FluxRiv[i][9]=0; 
    			} 
    		if(MD->Riv[i].LeftEle > 0)
    			{
			/*****************************************************************************/
			/* Lateral Surface Flux Calculation between River-Triangular element Follows */ 
			/*****************************************************************************/
			OLFeleToriv(MD->DummyY[MD->Riv[i].LeftEle - 1]+MD->Ele[MD->Riv[i].LeftEle - 1].zmax,MD->Ele[MD->Riv[i].LeftEle - 1].zmax,MD->Riv_Mat[MD->Riv[i].material-1].Cwr, MD->Riv[i].zmax,TotalY_Riv,MD->FluxRiv,i,2,MD->Riv[i].Length);
			/*********************************************************************************/
			/* Lateral Sub-surface Flux Calculation between River-Triangular element Follows */
			/*********************************************************************************/
			Dif_Y_Sub = (MD->DummyY[i+3*MD->NumEle] + MD->Riv[i].zmin) - (MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].LeftEle-1].zmin);
//			Avg_Y_Sub=(MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle-1].zmin-MD->Riv[i].zmin)>0?MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle-1].zmin-MD->Riv[i].zmin:0;
			/* This is head at river edge representation */
//			Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])-MD->Riv[i].zmin):0;
			/* This is head in neighboring cell represention */
			Avg_Y_Sub = MD->Ele[MD->Riv[i].LeftEle-1].zmin>MD->Riv[i].zmin?MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]:((MD->Ele[MD->Riv[i].LeftEle-1].zmin+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin?(MD->Ele[MD->Riv[i].LeftEle-1].zmin+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]-MD->Riv[i].zmin):0);
//			Avg_Y_Sub=avgY(MD->Riv[i].zmin,MD->Riv[i].zmin,MD->DummyY[i+3*MD->NumEle],Avg_Y_Sub);
			Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[i+3*MD->NumEle],Avg_Y_Sub);
			effK=MD->Riv[i].KsatH;
			Distance = sqrt(pow((MD->Riv[i].x - MD->Ele[MD->Riv[i].LeftEle - 1].x), 2) + pow((MD->Riv[i].y - MD->Ele[MD->Riv[i].LeftEle - 1].y), 2));                               
			Grad_Y_Sub = Dif_Y_Sub/Distance;
			/* take care of macropore effect */
			inabr=MD->Riv[i].LeftEle-1;
			AquiferDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
			effKnabr=effKH(MD->Ele[inabr].Macropore,MD->DummyY[inabr+2*MD->NumEle],AquiferDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
			Avg_Ksat=0.5*(effK+effKnabr);
			MD->FluxRiv[i][4]=MD->Riv[i].Length*Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub;
			/***********************************************************************************/
			/* Lateral Flux between rectangular element (beneath river) and triangular element */
			/***********************************************************************************/
                       	Dif_Y_Sub = (MD->DummyY[i+3*MD->NumEle+MD->NumRiv] + MD->Ele[i+MD->NumEle].zmin) - (MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].LeftEle-1].zmin);
//			Avg_Y_Sub=((MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle-1].zmin-MD->Riv[i].zmin)>0)?MD->Riv[i].zmin-MD->Ele[MD->Riv[i].LeftEle-1].zmin:MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle];
			/* This is head at river edge representation */
//			Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?MD->Riv[i].zmin-(MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)):MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle];
			/* This is head in neighboring cell represention */
			Avg_Y_Sub = MD->Ele[MD->Riv[i].LeftEle-1].zmin>MD->Riv[i].zmin?0:((MD->Ele[MD->Riv[i].LeftEle-1].zmin+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin?(MD->Riv[i].zmin-MD->Ele[MD->Riv[i].LeftEle-1].zmin):MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]);
//			Avg_Y_Sub=avgY(MD->Ele[i+MD->NumEle].zmin,MD->Ele[MD->Riv[i].LeftEle-1].zmin,MD->DummyY[i+3*MD->NumEle+MD->NumRiv],Avg_Y_Sub); 
			Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[i+3*MD->NumEle+MD->NumRiv],Avg_Y_Sub); 
			AquiferDepth=(MD->Ele[i+MD->NumEle].zmax-MD->Ele[i+MD->NumEle].zmin);                      
//			effK=MD->Ele[i+MD->NumEle].KsatH;
                     	effK=0.5*(effKH(MD->Ele[MD->Riv[i].LeftEle-1].Macropore,MD->DummyY[MD->Riv[i].LeftEle-1+2*MD->NumEle],MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin,MD->Ele[MD->Riv[i].LeftEle-1].macD,MD->Ele[MD->Riv[i].LeftEle-1].macKsatH,MD->Ele[MD->Riv[i].LeftEle-1].vAreaF,MD->Ele[MD->Riv[i].LeftEle-1].KsatH)+effKH(MD->Ele[MD->Riv[i].RightEle-1].Macropore,MD->DummyY[MD->Riv[i].RightEle-1+2*MD->NumEle],MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin,MD->Ele[MD->Riv[i].RightEle-1].macD,MD->Ele[MD->Riv[i].RightEle-1].macKsatH,MD->Ele[MD->Riv[i].RightEle-1].vAreaF,MD->Ele[MD->Riv[i].RightEle-1].KsatH));
                        inabr=MD->Riv[i].LeftEle-1;
			nabrAqDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
			effKnabr=effKH(MD->Ele[inabr].Macropore,MD->DummyY[inabr+2*MD->NumEle],nabrAqDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
			Avg_Ksat=0.5*(effK+effKnabr);
			Grad_Y_Sub = Dif_Y_Sub/Distance;                        /* take care of macropore effect */
			MD->FluxRiv[i][7]=MD->Riv[i].Length*Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub;
      			/* replace flux term */
      			for(j=0; j < 3; j++)
      				{
        			if(MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j] == MD->Riv[i].RightEle)
        				{
          				MD->FluxSurf[MD->Riv[i].LeftEle - 1][j] = -MD->FluxRiv[i][2];
					MD->FluxSub[MD->Riv[i].LeftEle - 1][j] =  -MD->FluxRiv[i][4];
					MD->FluxSub[MD->Riv[i].LeftEle - 1][j] = MD->FluxSub[MD->Riv[i].LeftEle - 1][j] -MD->FluxRiv[i][7];
          				break;
        				}
      				}      
    			}
    		if (MD->Riv[i].RightEle > 0)
    			{
			/*****************************************************************************/
			/* Lateral Surface Flux Calculation between River-Triangular element Follows */ 
			/*****************************************************************************/
			OLFeleToriv(MD->DummyY[MD->Riv[i].RightEle - 1]+MD->Ele[MD->Riv[i].RightEle - 1].zmax,MD->Ele[MD->Riv[i].RightEle - 1].zmax,MD->Riv_Mat[MD->Riv[i].material-1].Cwr, MD->Riv[i].zmax,TotalY_Riv,MD->FluxRiv,i,3,MD->Riv[i].Length);
			/*********************************************************************************/
                        /* Lateral Sub-surface Flux Calculation between River-Triangular element Follows */
			/*********************************************************************************/
                        Dif_Y_Sub = (MD->DummyY[i+3*MD->NumEle] + MD->Riv[i].zmin) - (MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].RightEle-1].zmin);
//			Avg_Y_Sub=(MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle-1].zmin-MD->Riv[i].zmin>0)?MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle-1].zmin-MD->Riv[i].zmin:0;
			/* This is head at river edge representation */
//			Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])-MD->Riv[i].zmin):0;
                        /* This is head in neighboring cell represention */
                        Avg_Y_Sub = MD->Ele[MD->Riv[i].RightEle-1].zmin>MD->Riv[i].zmin?MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]:((MD->Ele[MD->Riv[i].RightEle-1].zmin+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin?(MD->Ele[MD->Riv[i].RightEle-1].zmin+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]-MD->Riv[i].zmin):0);
//			Avg_Y_Sub=avgY(MD->Riv[i].zmin,MD->Riv[i].zmin,MD->DummyY[i+3*MD->NumEle],Avg_Y_Sub);
			Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[i+3*MD->NumEle],Avg_Y_Sub);
                        effK=MD->Riv[i].KsatH;
                        Distance = sqrt(pow((MD->Riv[i].x - MD->Ele[MD->Riv[i].RightEle - 1].x), 2) + pow((MD->Riv[i].y - MD->Ele[MD->Riv[i].RightEle - 1].y), 2));
                        Grad_Y_Sub = Dif_Y_Sub/Distance;
                        /* take care of macropore effect */
                        inabr=MD->Riv[i].RightEle-1;
                        AquiferDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
			effKnabr=effKH(MD->Ele[inabr].Macropore,MD->DummyY[inabr+2*MD->NumEle],AquiferDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
                        Avg_Ksat=0.5*(effK+effKnabr);
                        MD->FluxRiv[i][5]=MD->Riv[i].Length*Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub;
			/***********************************************************************************/
			/* Lateral Flux between rectangular element (beneath river) and triangular element */
			/***********************************************************************************/
                       	Dif_Y_Sub = (MD->DummyY[i+3*MD->NumEle+MD->NumRiv] + MD->Ele[i+MD->NumEle].zmin) - (MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].RightEle-1].zmin);
//			Avg_Y_Sub=((MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle-1].zmin-MD->Riv[i].zmin)>0)?MD->Riv[i].zmin-MD->Ele[MD->Riv[i].RightEle-1].zmin:MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle];
			/* This is head at river edge representation */
//			Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?MD->Riv[i].zmin-(MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)):MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle];
                        /* This is head in neighboring cell represention */                        
			Avg_Y_Sub = MD->Ele[MD->Riv[i].RightEle-1].zmin>MD->Riv[i].zmin?0:((MD->Ele[MD->Riv[i].RightEle-1].zmin+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin?(MD->Riv[i].zmin-MD->Ele[MD->Riv[i].RightEle-1].zmin):MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]);			
//			Avg_Y_Sub=avgY(MD->Ele[i+MD->NumEle].zmin,MD->Ele[MD->Riv[i].RightEle-1].zmin,MD->DummyY[i+3*MD->NumEle+MD->NumRiv],Avg_Y_Sub); 
			Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[i+3*MD->NumEle+MD->NumRiv],Avg_Y_Sub); 
			AquiferDepth=(MD->Ele[i+MD->NumEle].zmax-MD->Ele[i+MD->NumEle].zmin);                      
//			effK=MD->Ele[i+MD->NumEle].KsatH;
                     	effK=0.5*(effKH(MD->Ele[MD->Riv[i].LeftEle-1].Macropore,MD->DummyY[MD->Riv[i].LeftEle-1+2*MD->NumEle],MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin,MD->Ele[MD->Riv[i].LeftEle-1].macD,MD->Ele[MD->Riv[i].LeftEle-1].macKsatH,MD->Ele[MD->Riv[i].LeftEle-1].vAreaF,MD->Ele[MD->Riv[i].LeftEle-1].KsatH)+effKH(MD->Ele[MD->Riv[i].RightEle-1].Macropore,MD->DummyY[MD->Riv[i].RightEle-1+2*MD->NumEle],MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin,MD->Ele[MD->Riv[i].RightEle-1].macD,MD->Ele[MD->Riv[i].RightEle-1].macKsatH,MD->Ele[MD->Riv[i].RightEle-1].vAreaF,MD->Ele[MD->Riv[i].RightEle-1].KsatH));
                        inabr=MD->Riv[i].RightEle-1;
			nabrAqDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
			effKnabr=effKH(MD->Ele[inabr].Macropore,MD->DummyY[inabr+2*MD->NumEle],nabrAqDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
			Avg_Ksat=0.5*(effK+effKnabr);
			Grad_Y_Sub = Dif_Y_Sub/Distance;                        /* take care of macropore effect */
			MD->FluxRiv[i][8]=MD->Riv[i].Length*Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub;
      			/* replace flux item */
      			for(j=0; j < 3; j++)
      				{
        			if(MD->Ele[MD->Riv[i].RightEle - 1].nabr[j] == MD->Riv[i].LeftEle)
        				{
          				MD->FluxSurf[MD->Riv[i].RightEle - 1][j] = -MD->FluxRiv[i][3];
					MD->FluxSub[MD->Riv[i].RightEle - 1][j] =  -MD->FluxRiv[i][5];
					MD->FluxSub[MD->Riv[i].RightEle - 1][j] = MD->FluxSub[MD->Riv[i].RightEle - 1][j] -MD->FluxRiv[i][8];
          				break;
        				}
      				}
    			}
		Avg_Wid=CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,3);
		Dif_Y_Riv=(MD->Riv[i].zmin-(MD->DummyY[i + 3*MD->NumEle+MD->NumRiv]+MD->Ele[i+MD->NumEle].zmin))>0?MD->DummyY[i + 3*MD->NumEle]:MD->DummyY[i + 3*MD->NumEle]+MD->Riv[i].zmin-(MD->DummyY[i + 3*MD->NumEle+MD->NumRiv]+MD->Ele[i+MD->NumEle].zmin);
		Grad_Y_Riv=Dif_Y_Riv/MD->Riv[i].bedThick;
		MD->FluxRiv[i][6]=MD->Riv[i].KsatV*Avg_Wid*MD->Riv[i].Length*Grad_Y_Riv;
  		} 
		DEBUG FILE *tempfile;
		DEBUG tempfile=fopen("out.txt","ab");
	for(i=0; i<MD->NumEle; i++)
    		{   
     		for(j=0; j<3; j++)
      			{
        		DY[i] =  DY[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
        		DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle] - MD->FluxSub[i][j]/MD->Ele[i].area;
      			} 
      		DY[i+MD->NumEle] = DY[i+MD->NumEle]/(MD->Ele[i].Porosity*UNIT_C);	
      		DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle]/(MD->Ele[i].Porosity*UNIT_C);
		
		DEBUG //fprintf(tempfile,"A %lf\n", DY[i+0*MD->NumEle]);	
		DEBUG //fprintf(tempfile,"B %lf\n", DY[i+1*MD->NumEle]);	
		DEBUG //fprintf(tempfile,"C %lf\n", DY[i+2*MD->NumEle]);
		DY[i]=DY[i]/(UNIT_C);	
    		}
		DEBUG fclose(tempfile);  
   	for(i=0; i<MD->NumRiv; i++)
    		{
		for(j=0;j<=6;j++)
			{
			/* Note the limitation due to d(v)/dt=a*dy/dt+y*da/dt for CS other than rectangle */
			DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle]-MD->FluxRiv[i][j]/(MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3));
			}
		DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle]/(UNIT_C);
		DY[i+3*MD->NumEle+MD->NumRiv] = DY[i+3*MD->NumEle+MD->NumRiv] -MD->FluxRiv[i][7] -MD->FluxRiv[i][8]-MD->FluxRiv[i][9] -MD->FluxRiv[i][10]+MD->FluxRiv[i][6];
      		DY[i+3*MD->NumEle+MD->NumRiv] = DY[i+3*MD->NumEle+MD->NumRiv]/(MD->Ele[i+MD->NumEle].Porosity*MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3)*UNIT_C);
    		}






/* TRANSPORT */


/************************************************************/
/* IMPORTING THE HYDROLOGIC FLUXES FROM SOLVER COMPUTATIONS */
/************************************************************/

/* INITIALIZE HYDROLOGIC FLUXES: START */

/* ELEMENTS: START */
for (i=0; i<MD->NumEle; i++){
	for(j=0; j<3; j++){
		if(MD->Ele[i].nabr[j]>0){
			Dif_Y_Sub = (MD->DummyY[i+2*MD->NumEle] + MD->Ele[i].zmin) - (MD->DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle] + MD->Ele[MD->Ele[i].nabr[j]-1].zmin);
			Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[i+2*MD->NumEle],MD->DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle]);
			uSub[i][j]=Avg_Y_Sub*MD->Ele[i].edge[j]!=0?MD->FluxSub[i][j]/(Avg_Y_Sub*MD->Ele[i].edge[j]):0.0;
			uSub[i][j]=uSub[i][j]/UNIT_C;
			

			Dif_Y_Surf = (MD->SurfMode==1)?(MD->Ele[i].zmax-MD->Ele[MD->Ele[i].nabr[j] - 1].zmax):(MD->DummyY[i] + MD->Ele[i].zmax) - (MD->DummyY[MD->Ele[i].nabr[j] - 1] + MD->Ele[MD->Ele[i].nabr[j] - 1].zmax);
			Avg_Y_Surf=avgY(Dif_Y_Surf,MD->DummyY[i],MD->DummyY[MD->Ele[i].nabr[j]-1]);
			uSurf[i][j]=Avg_Y_Surf*MD->Ele[i].edge[j]!=0?MD->FluxSurf[i][j]/(Avg_Y_Surf*MD->Ele[i].edge[j]):0.0;
			uSurf[i][j]=uSurf[i][j]/UNIT_C;
		}
		else{
			uSub[i][j]=0.0;
			uSurf[i][j]=0.0;
		}
	}

	uPrep[i] = MD->EleNetPrep[i]/UNIT_C;
	if(MD->DummyY[i+2*MD->NumEle]>AquiferDepth-MD->Ele[i].infD){
		uInfilSat[i]=MD->EleViR[i]/UNIT_C;
		uInfilUsat[i]=0.0;
		MD->Recharge[i]=0.0;
		
		// *uETSurf, *uETUsat, *uETSat
		uETUsat[i]= 0.0;
		uETSat[i] = (MD->DummyY[i]<EPS/100)?MD->EleET[i][2]/UNIT_C:0.0;
	}
	else{
		uInfilSat[i]=0.0;
		uInfilUsat[i]=MD->EleViR[i]/UNIT_C;

		//
		uETUsat[i]= ((MD->DummyY[i]<EPS/100)&&(MD->DummyY[i+MD->NumEle]>EPS/100))?MD->EleET[i][2]/UNIT_C:0.0;
		uETSat[i] = 0.0;
	}
	uETSurf[i] = (MD->DummyY[i]<EPS/100)?0.0:MD->EleET[i][2]/UNIT_C;
}
for (i=0; i<MD->NumRiv; i++){
	for (j=0; j<3; j++){
		if(MD->Ele[MD->Riv[i].LeftEle-1].nabr[j]==MD->Riv[i].RightEle){
			uSub[MD->Riv[i].LeftEle-1][j]=0.0;
			uSurf[MD->Riv[i].LeftEle-1][j]=0.0;
		}
		if(MD->Ele[MD->Riv[i].RightEle-1].nabr[j]==MD->Riv[i].LeftEle){
			uSub[MD->Riv[i].RightEle-1][j]=0.0;
			uSurf[MD->Riv[i].RightEle-1][j]=0.0;
		}
	}
}
/* ELEMENTS: END */


/* STREAMS: START */
	// initialize qRiv, uRiv
	for(i=0; i<MD->NumRiv; i++){
		//for(j=0; j<=10; j++){
		for(j=0; j<11; j++){ //Yi - bug fix
			qRiv[i][j] = MD->FluxRiv[i][j]/UNIT_C;
		//	if(j==1) printf("%i %lf\t",i, qRiv[i][j]);
		}//if(i==1)printf("\n");

		if(MD->Riv[i].down > 0){
			CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,1);
			CrossAdown = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd,MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle],MD->Riv[MD->Riv[i].down - 1].coeff,1);
			As=0.5*(CrossA+CrossAdown); 
		
                      //  printf("%i\t",i);
                       // printf("%i %lf\t",i, qRiv[i][1]);getchar();
			uRiv[i][1] = As>0?qRiv[i][1]/As:0.0;

			uRiv[MD->Riv[i].down-1][0] += uRiv[i][1];
			uRivCount[i]++;			
		}
		else{
			switch(MD->Riv[i].down){
			case -3:
				As = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,1);
				uRiv[i][1] = As>0?qRiv[i][1]/As:0.0;
				break;
			default:
				printf("STREAM BOUNDARY CONDITION UNKNOWN\n");
				getchar();
				break;
			}
		}
	}
	for(i=0; i<MD->NumRiv; i++){
		if(uRivCount[i]>1)
			uRiv[i][0]=uRiv[i][0]/uRivCount[i];
	}
/* STREAMS: END */


/* INITIALIZE HYDROLOGIC FLUXES: END */



/************************************************************/
/* RHS CALCULATION FOR THE TRANSPORT COMPONENT OF THE MODEL */
/************************************************************/


/**********************************************************************************/
/* SURFACE     COMPONENT ONLY : START */

// * SURFACE-ADVECTION : START * /
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumEle; j++){
			for(k=0; k<3; k++){
				if(MD->Ele[j].nabr[k] > 0){
					if(uSurf[j][k]>0){
						Cs = DummyC[i*3*MD->NumEle + j];
						AgeCs = DummyAgeC[i*3*MD->NumEle + j];
					}
					else{
						Cs = DummyC[i*3*MD->NumEle + MD->Ele[j].nabr[k]-1];
						AgeCs = DummyAgeC[i*3*MD->NumEle + MD->Ele[j].nabr[k]-1];
					}
					Dif_Y_Surf = (MD->SurfMode==1)?(MD->Ele[j].zmax-MD->Ele[MD->Ele[j].nabr[k] - 1].zmax):(MD->DummyY[j] + MD->Ele[j].zmax) - (MD->DummyY[MD->Ele[j].nabr[k] - 1] + MD->Ele[MD->Ele[j].nabr[k] - 1].zmax);
					As=avgY(Dif_Y_Surf,MD->DummyY[j],MD->DummyY[MD->Ele[j].nabr[k]-1]);
					As = As * MD->Ele[j].edge[k];
					Ab = MD->Ele[j].area;
					if(MD->DummyY[j] > 0){
						DummyDC[i*3*MD->NumEle + j] += ((-uSurf[j][k]*As)/(-uSurf[j][k]*As+MD->DummyY[j]*Ab))*(Cs - DummyC[i*3*MD->NumEle + j]);
			if(fabs(DummyDC[i*3*MD->NumEle + j ]) > 5*24*60){ printf("E$1 %d",j); getchar;}
						AGE DummyAgeDC[i*3*MD->NumEle + j] += ((-uSurf[j][k]*As)/(-uSurf[j][k]*As+MD->DummyY[j]*Ab))*(AgeCs - DummyAgeC[i*3*MD->NumEle + j]);
					}
					else{
						DummyDC[i*3*MD->NumEle + j] += 0.0;
						DummyAgeDC[i*3*MD->NumEle + j] += 0.0;
					}
					if(isnan(DummyDC[i*3*MD->NumEle + j]) == 1){ printf("debug 1 %d %d %lf %lf %lf %lf %lf\n", i, j, DummyDC[i*3*MD->NumEle+j], uSurf[j][k], As, MD->DummyY[j], Ab);} 
				}
			}
		}
	}
// * SURFACE ADVECTION : END * /

// * SURFACE DISPERSION : START * /
	for(i=0; i<MD->NumEle; i++){
		Ux = 0.0;
		Uy = 0.0;
		uCount = 0;
		for(j=0; j<3; j++){
			if(MD->Ele[i].nabr[j] > 0){
				Ux += uSurf[i][j]*cos(SLOPE(MD->Ele[i].x, MD->Ele[i].y, MD->Ele[MD->Ele[i].nabr[j]-1].x, MD->Ele[MD->Ele[i].nabr[j]-1].y));
				Uy += uSurf[i][j]*sin(SLOPE(MD->Ele[i].x, MD->Ele[i].y, MD->Ele[MD->Ele[i].nabr[j]-1].x, MD->Ele[MD->Ele[i].nabr[j]-1].y));
				uCount++;
			}
		}
		Ux = Ux / uCount;
		Uy = Uy / uCount;
		U = sqrt(Ux*Ux + Uy*Uy);
		
		if(MD->AlphaInitModeSurf == 0)
			AlphaL = MD->AlphaLSurf; //??1.00;
		else{
			switch(MD->AlphaInitModeSurf){
			case 1:
				break;
			}
		}
				
                AlphaT = AlphaL/MD->AlphaTSurfFrac;
                Diffusivity = MD->DiffusionSurf;

		Disp[i][0] = U>0?AlphaL*Ux*Ux/U + AlphaT*Uy*Uy/U + 0.0:0.0;
                Disp[i][1] = U>0?AlphaL*Uy*Uy/U + AlphaT*Ux*Ux/U + 0.0:0.0;
		Disp[i][2] = U>0?(AlphaL - AlphaT) * Ux*Uy/U:0.0;
        }

        for(i=0; i<MD->NumSolute; i++){
                for(j=0; j<MD->NumEle; j++){
                        for(k=0; k<3; k++){
				if(MD->Ele[j].nabr[k]>0){
                                        Dxx = (Disp[j][0]+Disp[MD->Ele[j].nabr[k]-1][0])/2;
					Dyy = (Disp[j][1]+Disp[MD->Ele[j].nabr[k]-1][1])/2;
                                        Dxy = (Disp[j][2]+Disp[MD->Ele[j].nabr[k]-1][2])/2;
                                        Dyx = Dxy;
					// *if(j==2714 && MD->Ele[j].nabr[k]==3160){
                                        //        printf("Dxx=%lf Dxy=%lf Dyx=%lf Dyy=%lf\n", Dxx, Dxy, Dyx, Dyy); getchar();
                                        //}* /
                                        Cos = cos(SLOPE(MD->Ele[j].x, MD->Ele[j].y, MD->Ele[MD->Ele[j].nabr[k]-1].x, MD->Ele[MD->Ele[j].nabr[k]-1].y));
					Sin = sin(SLOPE(MD->Ele[j].x, MD->Ele[j].y, MD->Ele[MD->Ele[j].nabr[k]-1].x, MD->Ele[MD->Ele[j].nabr[k]-1].y));
					Dxx = Cos*Dxx + Sin*Dyx + Diffusivity; //Dxx=Diffusivity; //??
					Dxy = Cos*Dxy + Sin*Dyy;
					// *if(j==2714 && MD->Ele[j].nabr[k]==3160){
					//	printf("1Dxx=%lf Dxy=%lf Cos=%lf Sin=%lf\n", Dxx, Dxy, Cos, Sin); getchar();
					//}
					//if(j==215 && MD->Ele[j].nabr[k]==2708){
                                        //	printf("2Dxx=%lf Dxy=%lf Cos=%lf Sin=%lf\n", Dxx, Dxy, Cos, Sin); getchar(); 
                                        //}* /
                                        Dist = sqrt(pow((MD->Ele[j].x-MD->Ele[MD->Ele[j].nabr[k]-1].x),2)+pow((MD->Ele[j].y-MD->Ele[MD->Ele[j].nabr[k]-1].y),2));
                                        Grad_C = (DummyC[i*3*MD->NumEle + 0*MD->NumEle+j] - DummyC[i*3*MD->NumEle +0*MD->NumEle+MD->Ele[j].nabr[k]-1])/Dist;
					Grad_AgeC = (DummyAgeC[i*3*MD->NumEle + 0*MD->NumEle+j] - DummyAgeC[i*3*MD->NumEle +0*MD->NumEle+MD->Ele[j].nabr[k]-1])/Dist;
					
					Dif_Y_Surf = (MD->SurfMode==1)?(MD->Ele[j].zmax-MD->Ele[MD->Ele[j].nabr[k] - 1].zmax):(MD->DummyY[j] + MD->Ele[j].zmax) - (MD->DummyY[MD->Ele[j].nabr[k] - 1] + MD->Ele[MD->Ele[j].nabr[k] - 1].zmax);
					As=avgY(Dif_Y_Surf,MD->DummyY[j],MD->DummyY[MD->Ele[j].nabr[k]-1]);
                                        As = As * MD->Ele[j].edge[k];
                                        As =  MD->Ele[j].edge[k]; //Yi
					if(MD->DummyY[j] > EPS/100 && MD->DummyY[MD->Ele[j].nabr[k]-1] > EPS/100){
					DummyDC[i*3*MD->NumEle + 0*MD->NumEle+j] -= (Dxx+Dxy)*Grad_C*As/(MD->Ele[j].area);
	if(fabs(DummyDC[i*3*MD->NumEle + 0*MD->NumEle+j ]) > 5*24*60){ printf("E$2 %d",j); getchar;}
					AGE DummyAgeDC[i*3*MD->NumEle + 0*MD->NumEle+j] -= (Dxx+Dxy)*Grad_AgeC*As/(MD->Ele[j].area);
					}
					else{
						DummyDC[i*3*MD->NumEle + 0*MD->NumEle+j] += 0.0;
						DummyAgeDC[i*3*MD->NumEle + 0*MD->NumEle+j] += 0.0;
					}
					if(isnan(DummyDC[i*3*MD->NumEle + j]) == 1){ printf("debug 2 %d %d %lf %lf %lf %lf %lf %lf\n", i, j, DummyDC[i*3*MD->NumEle+j], Dxx, Dyy, Grad_C, As, MD->Ele[j].area);}
				}
			}
//Yi Surface Runoff flushing away DOC from the soil
                                        if (MD->DummyY[j] > EPS/100){
					DummyDC[i*3*MD->NumEle + 0*MD->NumEle+j] += (1/MD->DummyY[j]) * (MD->kovld/UNIT_C) * (MD->Ele[j].S0/MD->Kd - DummyC[i*3*MD->NumEle + 0*MD->NumEle+j]);
                                        }
// Yi
		}
	}
// * SURFACE DISPERSION : END   * /

// * SURFACE     COMPONENT ONLY : END   * /






/********************************************************************************/






// * UNSATURATED ZONE : START * /
/*
//Yi
        for (i=0; i<MD->NumSolute; i++){
                for(j=0; j<MD->NumEle; j++){
                        AquiferDepth = MD->Ele[j].zmax - MD->Ele[j].zmin;
                        Deficit = AquiferDepth - MD->DummyY[2*MD->NumEle + j];
                      //  Kd[i][j] = MD->Kd; // 0.03 from http://www.epa.gov/athens/learn2model/part-two/onsite/retard.html
                        if (MD->DummyY[2*MD->NumEle + j] > 0){
                           R[i][j] = (Deficit > MD->DummyY[1*MD->NumEle + j])? 1 + MD->Kd * MD->Bd*Deficit /(MD->Ele[j].Porosity*MD->DummyY[1*MD->NumEle + j]): 1 + MD->Kd * MD->Bd* MD->DummyY[1*MD->NumEle + j] /(MD->Ele[j].Porosity*MD->DummyY[1*MD->NumEle + j]) ;
                       }
                        else{
                           R[i][j] =  1 + MD->Kd * MD->Bd*AquiferDepth / (MD->Ele[j].Porosity *MD->DummyY[1*MD->NumEle + j]);
                       }
                }
        }
//Yi
*/
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumEle; j++){
			if(MD->Recharge[j] < 0){
				Cs = DummyC[i*3*MD->NumEle + 2*MD->NumEle + j];
				AgeCs = DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + j];
			}
			else{
				Cs = DummyC[i*3*MD->NumEle + 1*MD->NumEle + j];
				AgeCs = DummyAgeC[i*3*MD->NumEle + 1*MD->NumEle + j];
			}

			DummyDC[i*3*MD->NumEle + 1*MD->NumEle + j] += MD->DummyY[1*MD->NumEle + j]>0?((-MD->Recharge[j])/(-MD->Recharge[j]+MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*(Cs - DummyC[i*3*MD->NumEle + 1*MD->NumEle + j]):0.0;
//Yi
			if (j == cid){
                             tttt =  MD->DummyY[1*MD->NumEle + cid]>0?((-MD->Recharge[cid])/(-MD->Recharge[cid]+MD->DummyY[1*MD->NumEle + cid]*MD->Ele[cid].Porosity))*(Cs - DummyC[0*3*MD->NumEle + 1*MD->NumEle + cid]):0.0;
                             fprintf(omf[4],"%lf\t%lf\t",t,tttt) ;
                             tttt = Cs;
                        }
                        AquiferDepth=(MD->Ele[i].zmax-MD->Ele[i].zmin);
                        Deficit = AquiferDepth - MD->gwtb_old[j];

                        if (MD->DummyY[2*MD->NumEle + j]>0) {
                        DummyDC[i*3*MD->NumEle + 1*MD->NumEle + j] += MD->unsatY_old[j]> Deficit?(1/(MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*( MD->Bd*(MD->Tau/UNIT_C)*(DummyS[0*MD->NumEle + j]- MD->Kd*DummyC[1*MD->NumEle + j])* MD->unsatY_old[j]): (1/(MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*( MD->Bd*(MD->Tau/UNIT_C)*(DummyS[0*MD->NumEle + j]- MD->Kd*DummyC[1*MD->NumEle + j])*Deficit) ;
                           if (j ==cid){
                        tttt1 = MD->unsatY_old[j]> Deficit?(1/(MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*( MD->Bd*(MD->Tau/UNIT_C)*(DummyS[0*MD->NumEle + j]- MD->Kd*DummyC[1*MD->NumEle + j])* MD->unsatY_old[j]): (1/(MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*( MD->Bd*(MD->Tau/UNIT_C)*(DummyS[0*MD->NumEle + j]- MD->Kd*DummyC[1*MD->NumEle + j])*Deficit) ;
                        fprintf(omf[4],"%lf\t",tttt1); 
                            }
                        }else{
                        DummyDC[i*3*MD->NumEle + 1*MD->NumEle + j] += (1/(MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*( MD->Bd*(MD->Tau/UNIT_C)*(DummyS[0*MD->NumEle + j]- MD->Kd*DummyC[1*MD->NumEle + j])* AquiferDepth);
                           if (j ==cid){
                        tttt1 = (1/(MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*( MD->Bd*(MD->Tau/UNIT_C)*(DummyS[0*MD->NumEle + j]- MD->Kd*DummyC[1*MD->NumEle + j])* AquiferDepth);
                        fprintf(omf[4],"%lf\t",tttt1);
                            }
                        }
//printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",j,tttt1,DummyS[0*MD->NumEle + j],DummyC[1*MD->NumEle + j],MD->unsatY_old[j],Deficit);
                        
// Surface Runoff flushing away DOC from the soil
                                        if (MD->DummyY[j] > EPS/100){
					   DummyDC[i*3*MD->NumEle + 1*MD->NumEle+j] += -(1/(MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity)) * (MD->kovld/UNIT_C) * (MD->Ele[j].S0/MD->Kd - DummyC[i*3*MD->NumEle + 0*MD->NumEle+j]);
                                           if (j==cid){
                                              fprintf(omf[4],"%lf\t",  -(1/(MD->DummyY[1*MD->NumEle + cid]*MD->Ele[cid].Porosity)) * (MD->kovld/UNIT_C) * (MD->Ele[cid].S0/MD->Kd - DummyC[i*3*MD->NumEle + 0*MD->NumEle+cid]));
                                           }
                                        }else{
                                           if (j==cid){
                                              fprintf(omf[4],"%lf\t", 0.0);
                                           }
                                        }
// Yi
                     
			DummyDC[i*3*MD->NumEle + 2*MD->NumEle + j] += MD->DummyY[2*MD->NumEle + j]>0?(MD->Recharge[j]/(MD->Recharge[j]+MD->DummyY[2*MD->NumEle + j]*MD->Ele[j].Porosity))*(Cs - DummyC[i*3*MD->NumEle + 2*MD->NumEle + j]):0.0;

			AGE{
				DummyAgeDC[i*3*MD->NumEle + 1*MD->NumEle + j] += MD->DummyY[1*MD->NumEle + j]>0?((-MD->Recharge[j])/(-MD->Recharge[j]+MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity))*(AgeCs - DummyAgeC[i*3*MD->NumEle + 1*MD->NumEle + j]):0.0;
				DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle + j] += MD->DummyY[2*MD->NumEle + j]>0?(MD->Recharge[j]/(MD->Recharge[j]+MD->DummyY[2*MD->NumEle + j]*MD->Ele[j].Porosity))*(AgeCs - DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + j]):0.0;
			}
		}
	}
                         fprintf(omf[4],"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",MD->Recharge[cid],tttt,MD->Recharge[cid]*Cs*MD->Ele[cid].area, DummyC[0*3*MD->NumEle + 1*MD->NumEle + cid], DummyDC[0*3*MD->NumEle + 1*MD->NumEle + cid], DummyDC[0*3*MD->NumEle + 2*MD->NumEle + cid]); //Yi
                         fprintf(omf[4],"%lf\t%lf\t",MD->DummyY[1*MD->NumEle + cid],MD->DummyY[2*MD->NumEle + cid]); //Yi

// Yi
   
// * UNSATURATED ZONE : END   * /






/********************************************************************************/

// * RETARDATION : START * //

realtype bulkdensity;
//bulkdensity = 1.99; // TODO this can be improved by adding this as an input in the soil file

	for (i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumEle; j++){
			// LINEAR SORPTION ISOTHERM
//			Kd[i][j] = 0.05; // 0.03 from http://www.epa.gov/athens/learn2model/part-two/onsite/retard.html
			// FREUNDLICH SORPTION ISOTHERM

			// LANGMUIR SORPTION ISOTHERM

//			R[i][j] = 1 + Kd[i][j] * (bulkdensity / MD->Ele[j].Porosity);
//			R[i][j] = 2.9;
                        //Yi
                        R[i][j] = 1.0 ; 
                        //Yi
		}
	}

// * RETARDATION : END * //



/* SUB-SURFACE COMPONENT ONLY : START */

/* ADVECTION : START */
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumEle; j++){
			for(k=0; k<3; k++){
				if(MD->Ele[j].nabr[k] > 0){
					//if(fabs(uSub[j][k])>0) printf("Advection: %lf %d %d %lf\n", t, j, k,uSub[j][k]);
					if(uSub[j][k]>0){
						Cs = DummyC[i*3*MD->NumEle + 2*MD->NumEle+j];
						AgeCs = DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle+j];
					}						
					else{
						Cs = DummyC[i*3*MD->NumEle + 2*MD->NumEle+MD->Ele[j].nabr[k]-1];
						AgeCs = DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle+MD->Ele[j].nabr[k]-1];
					}						
					//if(Cs == 1.0)
                                          //      printf("u=%lf, i=%d, j=%d, k=%d\n", uSub[j][k], i, j, k);					
	
					Dif_Y_Sub = (MD->DummyY[j+2*MD->NumEle] + MD->Ele[i].zmin) - (MD->DummyY[MD->Ele[j].nabr[k]-1 + 2*MD->NumEle] + MD->Ele[MD->Ele[j].nabr[k]-1].zmin);
					Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[j+2*MD->NumEle],MD->DummyY[MD->Ele[j].nabr[k]-1 + 2*MD->NumEle]);
					As = Avg_Y_Sub * MD->Ele[j].edge[k];
					Ab = MD->Ele[j].area;
					if(MD->DummyY[2*MD->NumEle+j] > 0){
						DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j] += (1/R[i][j])*((-uSub[j][k]*As)/(-uSub[j][k]*As+MD->DummyY[2*MD->NumEle+j]*Ab*MD->Ele[j].Porosity))*(Cs - DummyC[i*3*MD->NumEle + 2*MD->NumEle+j]);
	if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j ]) > 5*24*60){ printf("E$3 %d",j); getchar;}
						AGE DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle+j] += (1/R[i][j])*((-uSub[j][k]*As)/(-uSub[j][k]*As+MD->DummyY[2*MD->NumEle+j]*Ab*MD->Ele[j].Porosity))*(AgeCs - DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle+j]);
					}
					//if(j==9)
					//	{printf("a j=%d %lf %lf %lf %lf",j,t, Cs, DummyC[i*3*MD->NumEle + 2*MD->NumEle+j],DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j]); getchar();}					
if(isnan(DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j]) == 1){ printf("debug 3 %d %d %lf %lf %lf %lf\n", i, j, uSub[j][k], As, Ab, DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j]);}
					//if(j==304 && Cs == 1.0)
					//	printf("304->%lf\n", DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j]);
					//if(j==305 && Cs == 1.0)
					//	printf("305->%lf\n", DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j]);
				}
			}
		}
	}
// * ADVECTION : END * /

// * DISPERSION : START * /
//if(MD->NumSolute>0)
	for(i=0; i<MD->NumEle; i++){
		Ux = 0.0;
		Uy = 0.0;
		uCount = 0;Distance = 0;
		for(j=0; j<3; j++){
			if(MD->Ele[i].nabr[j] > 0){
				Ux += uSub[i][j]*cos(SLOPE(MD->Ele[i].x, MD->Ele[i].y, MD->Ele[MD->Ele[i].nabr[j]-1].x, MD->Ele[MD->Ele[i].nabr[j]-1].y));
				Uy += uSub[i][j]*sin(SLOPE(MD->Ele[i].x, MD->Ele[i].y, MD->Ele[MD->Ele[i].nabr[j]-1].x, MD->Ele[MD->Ele[i].nabr[j]-1].y));
				uCount++;
				Distance += sqrt(pow((MD->Ele[i].x - MD->Ele[MD->Ele[i].nabr[j] - 1].x), 2) + pow((MD->Ele[i].y - MD->Ele[MD->Ele[i].nabr[j] - 1].y), 2));
			}
		}
		Ux = Ux / uCount;
		Uy = Uy / uCount; Distance = Distance / uCount;
		U = sqrt(Ux*Ux + Uy*Uy);
		if(MD->AlphaInitModeSat==0)
			AlphaL=MD->AlphaLSat;
		else{
			switch(MD->AlphaInitModeSat){
			case 1:
				if(Distance>1.0)
					AlphaL = 0.83*pow(log10(Distance), 2.414);//??0.01;
				else
					AlphaL = 0.01;
				break;
			}
		}
		//??printf("\n Distance AlphaL = %lf %lf", Distance, AlphaL);
		AlphaT = AlphaL/MD->AlphaTSatFrac;
		Diffusivity = MD->DiffusionSat; //??0.01000;

		//Disp[i][0] = U>0?AlphaL*Ux*Ux/U + AlphaT*Uy*Uy/U + Diffusivity:Diffusivity;
		//Disp[i][1] = U>0?AlphaL*Uy*Uy/U + AlphaT*Ux*Ux/U + Diffusivity:Diffusivity;
		Disp[i][0] = U>0?AlphaL*Ux*Ux/U + AlphaT*Uy*Uy/U + 0.0:0.0;
                Disp[i][1] = U>0?AlphaL*Uy*Uy/U + AlphaT*Ux*Ux/U + 0.0:0.0;
		Disp[i][2] = U>0?(AlphaL - AlphaT) * Ux*Uy/U:0.0;
	}//??getchar();

	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumEle; j++){
                      
                     AquiferDepth=(MD->Ele[j].zmax-MD->Ele[j].zmin); //Yi
			for(k=0; k<3; k++){
				if(MD->Ele[j].nabr[k]>0){
					Dxx = (Disp[j][0]+Disp[MD->Ele[j].nabr[k]-1][0])/2;
					Dyy = (Disp[j][1]+Disp[MD->Ele[j].nabr[k]-1][1])/2;
					Dxy = (Disp[j][2]+Disp[MD->Ele[j].nabr[k]-1][2])/2;
					Dyx = Dxy;
					
					Cos = cos(SLOPE(MD->Ele[j].x, MD->Ele[j].y, MD->Ele[MD->Ele[j].nabr[k]-1].x, MD->Ele[MD->Ele[j].nabr[k]-1].y));
					Sin = sin(SLOPE(MD->Ele[j].x, MD->Ele[j].y, MD->Ele[MD->Ele[j].nabr[k]-1].x, MD->Ele[MD->Ele[j].nabr[k]-1].y));
					//if(Cos*Dxx + Sin*Dyx != 0.0){ printf("dispersion not zero\n"); getchar(); }	
					Dxx = Cos*Dxx + Sin*Dyx + Diffusivity; //??Dxx=Diffusivity;
					Dxy = Cos*Dxy + Sin*Dyy;
					//printf("Diff = %d %d %lf\n", j, k, Dxx);
					Dist = sqrt(pow((MD->Ele[j].x-MD->Ele[MD->Ele[j].nabr[k]-1].x),2)+pow((MD->Ele[j].y-MD->Ele[MD->Ele[j].nabr[k]-1].y),2));
					Grad_C = (DummyC[i*3*MD->NumEle + 2*MD->NumEle+j] - DummyC[i*3*MD->NumEle +2*MD->NumEle+MD->Ele[j].nabr[k]-1])/Dist;
					Grad_AgeC = (DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle+j] - DummyAgeC[i*3*MD->NumEle +2*MD->NumEle+MD->Ele[j].nabr[k]-1])/Dist;

					Dif_Y_Sub = (MD->DummyY[j+2*MD->NumEle] + MD->Ele[j].zmin) - (MD->DummyY[MD->Ele[j].nabr[k]-1 + 2*MD->NumEle] + MD->Ele[MD->Ele[j].nabr[k]-1].zmin);
					Avg_Y_Sub=avgY(Dif_Y_Sub,MD->DummyY[j+2*MD->NumEle],MD->DummyY[MD->Ele[j].nabr[k]-1 + 2*MD->NumEle]);
                                        As = Avg_Y_Sub * MD->Ele[j].edge[k];
                                        As = MD->Ele[j].edge[k];//Yi
					if(MD->DummyY[2*MD->NumEle+j] > 0 && MD->DummyY[2*MD->NumEle+MD->Ele[j].nabr[k]-1] > 0){
					DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j] -= (1/R[i][j])*(Dxx+Dxy)*Grad_C*As/MD->Ele[j].area;
		//			DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j] -= (1/R[i][j])*(Dxx+Dxy)*Grad_C*As/(MD->Ele[j].area*Avg_Y_Sub);
	if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j ]) > 5*24*60){ printf("E$4 %d,j"); getchar;}
					AGE DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle+j] -= (1/R[i][j])*(Dxx+Dxy)*Grad_AgeC*As/MD->Ele[j].area;
					}	
if(isnan(DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j]) == 1){ printf("debug 4 %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j]);}
					//if(fabs(Dxx*Grad_C*As/MD->Ele[j].area)>0.0){
						//printf("\nDiffusion: %d %d %.15lf %lf", j, MD->Ele[j].nabr[k], Dxx*Grad_C*As/MD->Ele[j].area, t);
					//	getchar();
					//}
				}
			}

//Yi
//   realtype test;
                   if (MD->DummyY[2*MD->NumEle + j] > 0.0000){
                       if (MD->gwtb_old[j] < AquiferDepth - MD->unsatY_old[j]){
                     //  if (MD->DummyY[2*MD->NumEle + j] < AquiferDepth){
                      DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j] -= (1.0/(MD->Ele[j].Porosity*MD->DummyY[2*MD->NumEle +j]))*(- MD->Bd*MD->gwtb_old[j]*(MD->Tau/UNIT_C)*(DummyS[1*MD->NumEle + j] - MD->Kd*DummyC[2*MD->NumEle + j])); 
                           if (j == cid){
                      tttt = -(1.0/(MD->Ele[cid].Porosity*MD->DummyY[2*MD->NumEle +cid]))*(- MD->Bd*MD->gwtb_old[cid]*(MD->Tau/UNIT_C)*(DummyS[1*MD->NumEle + cid] - MD->Kd*DummyC[2*MD->NumEle + cid]));
                      fprintf(omf[4],"%lf\t",tttt);
                           }
                        }
                      else{
                      DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j] -= (1.0/(MD->Ele[j].Porosity*MD->DummyY[2*MD->NumEle +j]))*(- MD->Bd*(AquiferDepth - MD->unsatY_old[j])*(MD->Tau/UNIT_C)*(DummyS[1*MD->NumEle + j] - MD->Kd*DummyC[2*MD->NumEle + j]));
		          if (j == cid){
                      tttt = -(1.0/(MD->Ele[cid].Porosity*MD->DummyY[2*MD->NumEle +cid]))*(- MD->Bd*(AquiferDepth - MD->unsatY_old[cid])*(MD->Tau/UNIT_C)*(DummyS[1*MD->NumEle + cid] - MD->Kd*DummyC[2*MD->NumEle + cid]));
                      fprintf(omf[4],"%lf\t",tttt);
                           }
 
                     // DummyDC[i*3*MD->NumEle + 2*MD->NumEle+j] -= (1.0/(MD->Ele[j].Porosity*MD->DummyY[2*MD->NumEle +j]))*(- MD->Bd*(AquiferDepth )*(MD->Tau/UNIT_C)*(DummyS[1*MD->NumEle + j] - MD->Kd*DummyC[2*MD->NumEle + j]));
                      }
                   }else{
		          if (j == cid){
                      tttt = 0.0;
                      fprintf(omf[4],"%lf\t",tttt);
                           }
                   }
//Yi

		}
                  fprintf(omf[4],"%lf\t%lf\t",- MD->Bd*MD->gwtb_old[cid]*(MD->Tau/UNIT_C),(DummyS[1*MD->NumEle + cid] - MD->Kd*DummyC[2*MD->NumEle + cid]));
	}
// * DISPERSION : END * /

// * SUB-SURFACE COMPONENT ONLY : END   * /






/**********************************************************************************/
/* RIVER     COMPONENT ONLY : START */



// * ADVECTION: START * /
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumRiv; j++){
			Ab =(MD->Riv[j].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,MD->DummyY[j + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,3));

			// * RIVER & RIVER + DOWN : START * /
			if(qRiv[j][1]>0){
				Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
				AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
			}				
			else{
				if(MD->Riv[j].down>0){
					Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down -1];
					AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down -1];
				}
				else{
					//BC ??
				}
			}
		
			DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][1])/(-qRiv[j][1]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j ]) > 5*24*60){ printf("E$5 %d",j); getchar;}

DEBUG
//tDC = MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][1])/(-qRiv[j][1]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
//if(tDC < 0) {printf("D1 %lf\n", tDC); getchar();}

			AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][1])/(-qRiv[j][1]+Ab*MD->DummyY[3*MD->NumEle+j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;

			if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("debug 5 %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]);}
			
			if(MD->Riv[j].down>0){
				Ab_down =(MD->Riv[MD->Riv[j].down-1].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].interpOrd,MD->DummyY[MD->Riv[j].down-1 + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].coeff,3));

				DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1] += MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]>0?((qRiv[j][1])/(qRiv[j][1]+Ab_down*MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1 ]) > 5*24*60){ printf("E$6 %d",j); getchar;}

DEBUG
//tDC = MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]>0?((qRiv[j][1])/(qRiv[j][1]+Ab_down*MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]):0.0;
//if(tDC < 0) {printf("D2 %lf\n", tDC); getchar();}

				AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1] += MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]>0?((qRiv[j][1])/(qRiv[j][1]+Ab_down*MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]):0.0;

				if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]) == 1){ printf("debug 6 %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]);}
			}
			// * RIVER & RIVER + DOWN : END * /










			// * RIVER & SURFACE : START * /
			if(qRiv[j][2] > 0){
				Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
				AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
			}				
			else{
				Cs = DummyC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle - 1];
				AgeCs = DummyAgeC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle - 1];
			}
			if(qRiv[j][2] < 0)				
			DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][2])/(-qRiv[j][2]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j ]) > 5*24*60){ printf("E$7 %d",j); getchar;}

DEBUG
//tDC =  MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][2])/(-qRiv[j][2]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
//if(tDC < 0) {printf("D3 %lf %lf %lf %lf %lf %lf %lf\n", tDC, -qRiv[j][2], Ab*MD->DummyY[3*MD->NumEle+j], Cs, DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j],DummyC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle - 1], MD->DummyY[MD->Riv[j].LeftEle - 1]); getchar();}

	if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("debug 7 %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]);}
			if(qRiv[j][2] > 0)
			DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1] += MD->DummyY[0*MD->NumEle + MD->Riv[j].LeftEle-1]>0?((qRiv[j][2])/(qRiv[j][2]+MD->Ele[MD->Riv[j].LeftEle-1].area*MD->DummyY[0*MD->NumEle + MD->Riv[j].LeftEle-1]))*(Cs-DummyC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1 ]) > 5*24*60){ printf("E$8 %d",j); getchar;}

DEBUG
//tDC = MD->DummyY[0*MD->NumEle + MD->Riv[j].LeftEle-1]>0?((qRiv[j][2])/(qRiv[j][2]+MD->Ele[MD->Riv[j].LeftEle-1].area*MD->DummyY[0*MD->NumEle + MD->Riv[j].LeftEle-1]))*(Cs-DummyC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1]):0.0;
//if(tDC < 0) {printf("D4 %lf\n", tDC); getchar();}

	if(isnan(DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1]) == 1){ printf("debug 8 %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1]);}
		AGE{
		    if(qRiv[j][2]<0)
			DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][2])/(-qRiv[j][2]+Ab*MD->DummyY[3*MD->NumEle+j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
		    if(qRiv[j][2]>0)
			DummyAgeDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1] += MD->DummyY[0*MD->NumEle + MD->Riv[j].LeftEle-1]>0?((qRiv[j][2])/(qRiv[j][2]+MD->Ele[MD->Riv[j].LeftEle-1].area*MD->DummyY[0*MD->NumEle + MD->Riv[j].LeftEle-1]))*(AgeCs-DummyAgeC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].LeftEle-1]):0.0;
		}

			if(qRiv[j][3] > 0){
                                Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
                                AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
			}				
                        else{
                                Cs = DummyC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle - 1];
                                AgeCs = DummyAgeC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle - 1];
			}				
			if(i==0 && (qRiv[j][2]>0 || qRiv[j][3]>0)) printf("%d %lf %lf\n", j, qRiv[j][2], qRiv[j][3]); 
                        DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][3])/(-qRiv[j][3]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j ]) > 5*24*60){ printf("E$9 %d",j); getchar;}

DEBUG
//tDC =  MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][3])/(-qRiv[j][3]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
//if(tDC < 0) {printf("D5 %lf\n", tDC); getchar();}

                        DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1] += MD->DummyY[0*MD->NumEle + MD->Riv[j].RightEle-1]>0?((qRiv[j][3])/(qRiv[j][3]+MD->Ele[MD->Riv[j].RightEle-1].area*MD->DummyY[0*MD->NumEle + MD->Riv[j].RightEle-1]))*(Cs-DummyC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1]):0.0;

DEBUG
//tDC = MD->DummyY[0*MD->NumEle + MD->Riv[j].RightEle-1]>0?((qRiv[j][3])/(qRiv[j][3]+MD->Ele[MD->Riv[j].RightEle-1].area*MD->DummyY[0*MD->NumEle + MD->Riv[j].RightEle-1]))*(Cs-DummyC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1]):0.0;
//if(tDC < 0) {printf("D6 %lf\n", tDC); getchar();}


if(fabs(DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1 ]) > 5*24*60){ printf("E$10 %d",j); getchar;}
			if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("debug 7a %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]);}
			if(isnan(DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1]) == 1){ printf("debug 8a %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1]);}

			AGE{
				 DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][3])/(-qRiv[j][3]+Ab*MD->DummyY[3*MD->NumEle+j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
				 DummyAgeDC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1] += MD->DummyY[0*MD->NumEle + MD->Riv[j].RightEle-1]>0?((qRiv[j][3])/(qRiv[j][3]+MD->Ele[MD->Riv[j].RightEle-1].area*MD->DummyY[0*MD->NumEle + MD->Riv[j].RightEle-1]))*(AgeCs-DummyAgeC[i*3*MD->NumEle + 0*MD->NumEle + MD->Riv[j].RightEle-1]):0.0;
			}
			// * RIVER & SURFACE : END   * /








			// * RIVER & SUB-SURFACE : START   * /
			if(qRiv[j][4] > 0){
                                Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
                                AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
			}				
                        else{
                                Cs = DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle - 1];
                                AgeCs = DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle - 1];
			}				
                        DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][4])/(-qRiv[j][4]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j ]) > 5*24*60){ printf("E$11 %d",j); getchar;}
                        DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]>0?((qRiv[j][4])/(qRiv[j][4]+MD->Ele[MD->Riv[j].LeftEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]*MD->Ele[MD->Riv[j].LeftEle-1].Porosity))*(Cs-DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1 ]) > 5*24*60){ printf("E$12 %d",j); getchar;}
			if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("debug 7b %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]);}                         
			if(isnan(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]) == 1){ printf("debug 8b %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]);}
		AGE{
			DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][4])/(-qRiv[j][4]+Ab*MD->DummyY[3*MD->NumEle+j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
                        DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]>0?((qRiv[j][4])/(qRiv[j][4]+MD->Ele[MD->Riv[j].LeftEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]*MD->Ele[MD->Riv[j].LeftEle-1].Porosity))*(AgeCs-DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]):0.0;
		}

                        if(qRiv[j][5] > 0){
                                Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
                                AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
			}				
                        else{
                                Cs = DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle - 1];
                                AgeCs = DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle - 1];
			}				
                        DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][5])/(-qRiv[j][5]+Ab*MD->DummyY[3*MD->NumEle+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) > 5*24*60){ printf("E$13 %d",j); getchar;}
                        DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]>0?((qRiv[j][5])/(qRiv[j][5]+MD->Ele[MD->Riv[j].RightEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]*MD->Ele[MD->Riv[j].RightEle-1].Porosity))*(Cs-DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]) > 5*24*60){ printf("E$14 %d",j); getchar;}
			if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("debug 7c %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]);}
			if(isnan(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]) == 1){ printf("debug 8c %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]);}
		AGE{
			DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+j]>0?((-qRiv[j][5])/(-qRiv[j][5]+Ab*MD->DummyY[3*MD->NumEle+j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;
                        DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]>0?((qRiv[j][5])/(qRiv[j][5]+MD->Ele[MD->Riv[j].RightEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]*MD->Ele[MD->Riv[j].RightEle-1].Porosity))*(AgeCs-DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]):0.0;
		}


			// * RIVER & SUB-SURFACE : END   * /







		}
	}
// * ADVECTION: END * /






// * DISPERSION: START */
/*	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumRiv; j++){
			Ux = (uRiv[j][0] + uRiv[j][1])/2;
			U  = sqrt(Ux*Ux);
			if(MD->AlphaInitModeRiv == 0)
				AlphaL=MD->AlphaLRiv;
			else{
				switch(MD->AlphaInitModeRiv){
				case 1:
					HydraulicRadius = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,MD->DummyY[j + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,2)/CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,MD->DummyY[j + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,2);
					uShear = sqrt(fabs(GRAV*HydraulicRadius*MD->Riv[j].BedSlope));            // * U* = sqrt(gRS) * /
					AlphaL = fabs(10.612*MD->DummyY[3*MD->NumEle+j]*(U/uShear));//??0.01;
					break;
				case 2:
					uShear = sqrt(fabs(GRAV*MD->DummyY[3*MD->NumEle+j]*MD->Riv[j].BedSlope)); // * U* = sqrt(gHS) * /
					AlphaL = fabs(10.612*MD->DummyY[3*MD->NumEle+j]*(U/uShear));//??0.01;
					break;
				}
			}
			//??if(t>43150) printf("AlphaRiv\t%lf\t%lf\t%lf\t%lf\n", DummyY[3*MD->NumEle+j], U, uShear, AlphaL);
			AlphaT = AlphaL/MD->AlphaTRivFrac;
			Diffusivity = MD->DiffusionRiv; //??0.00000;
			//??Ux = (uRiv[j][0] + uRiv[j][1])/2;
			//??U  = sqrt(Ux*Ux);
			DispRiv[j] = (Ux!=0 && uShear!=0)?AlphaL*Ux*Ux/U:0.0;
		}
	}//??getchar();
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumRiv; j++){
			if(MD->Riv[j].down>0){
				Dxx = (DispRiv[j] + DispRiv[MD->Riv[j].down-1])/2+Diffusivity;
				Grad_C = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] - DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down -1];
				AGE Grad_AgeC = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] - DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down -1];

				CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->DummyY[i + 3*MD->NumEle],MD->Riv[i].coeff,1);
				CrossAdown = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd,MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle],MD->Riv[MD->Riv[i].down - 1].coeff,1);
				As=0.5*(CrossA+CrossAdown);
				//As = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,Avg_Y_Riv,MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,1);
				Ab =(MD->Riv[j].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,MD->DummyY[j + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,3));
				if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("debug 9a %d %d %lf %lf %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j], Dxx, MD->DummyY[3*MD->NumEle+j]);}

				if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]) == 1){ printf("debug 10a %d %d %lf %lf %lf %lf %lf %lf - %lf %lf %lf %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1], Dxx, Grad_C, As, Ab, MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1], DispRiv[j], DispRiv[MD->Riv[j].down-1], uRiv[j][0], uRiv[j][1]);}

				if(MD->DummyY[3*MD->NumEle+j]>0 && MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]>0){
					DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] -= MD->DummyY[3*MD->NumEle+j]>0?Dxx*Grad_C*As/(Ab*MD->DummyY[3*MD->NumEle+j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) > 5*24*60){ printf("E$15 %d",j); getchar;}
					AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] -= MD->DummyY[3*MD->NumEle+j]>0?Dxx*Grad_AgeC*As/(Ab*MD->DummyY[3*MD->NumEle+j]):0.0;					
					//??As = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].interpOrd,Avg_Y_Riv,MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].coeff,1);
					Ab =(MD->Riv[MD->Riv[j].down-1].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].interpOrd,MD->DummyY[MD->Riv[j].down-1 + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].coeff,3));
					DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1] += MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]>0?Dxx*Grad_C*As/(Ab*MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]) > 5*24*60){ printf("E$16 %d",j); getchar;}
					AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1] += MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]>0?Dxx*Grad_AgeC*As/(Ab*MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1]):0.0;				
					if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("debug 9 %d %d %lf %lf %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*MD->NumRiv + j], Dxx, MD->DummyY[3*MD->NumEle+j]);}
					if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1]) == 1){ printf("debug 10a %d %d %lf %lf %lf %lf %lf %lf - %lf %lf %lf %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->Riv[j].down-1], Dxx, Grad_C, As, Ab, MD->DummyY[3*MD->NumEle+MD->Riv[j].down-1], DispRiv[j], DispRiv[MD->Riv[j].down-1], uRiv[j][0], uRiv[j][1]);}
				}
			}
			else{
				//TODO BC //??
			}
		}
	}*/


// * DISPERSION: END */



// * RIVER     COMPONENT ONLY : START * /







/**********************************************************************************/



// * PRECIPITATION + ET + INFIL : START * /
//??    CPrep = Interpolation(&MD->TSD_PChem[MD->Ele[i].prep-1], t); //1.0;
    for(i=0; i<MD->NumSolute; i++){
	for(j=0; j<MD->NumEle; j++){
		CPrep = Interpolation(&MD->TSD_PChem[i*MD->NumPrep+MD->Ele[j].prep-1], t);
		AgeCPrep = 0;
		//printf("CPrep= %lf\t", CPrep);
		CInfil = MD->DummyY[j]>0?DummyC[i*3*MD->NumEle + j]:CPrep;
		AgeCInfil= MD->DummyY[j]>0?DummyAgeC[i*3*MD->NumEle + j]:AgeCPrep;
		//CInfil = CPrep;
		//printf("CPrep= %lf\t", CInfil);

		if(isnan(DummyDC[i*3*MD->NumEle + j]) == 1){ printf("debug 11 %d %d %lf %lf %lf %lf\n", i, j, DummyDC[i*3*MD->NumEle + j], MD->DummyY[j], uInfilUsat[i], uInfilSat[i]);}
		DummyDC[i*3*MD->NumEle + j] += (uPrep[j]+MD->DummyY[j])>0.0?(uPrep[j]/(uPrep[j]+MD->DummyY[j])) * (CPrep  - DummyC[i*3*MD->NumEle + j]):0.0;
		AGE DummyAgeDC[i*3*MD->NumEle + j] += (uPrep[j]+MD->DummyY[j])>0.0?(uPrep[j]/(uPrep[j]+MD->DummyY[j])) * (AgeCPrep  - DummyAgeC[i*3*MD->NumEle + j]):0.0;
		//DummyDC[i*3*MD->NumEle + j] += (MD->DummyY[j])>0.0?(DummyC[i*3*MD->NumEle + j]<=0?CPrep:(uPrep[j]/(uPrep[j]+MD->DummyY[j])) * (CPrep  - DummyC[i*3*MD->NumEle + j])):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + j]) > 5*24*60){ printf("E$17 %d",j); getchar;}
//		if(j==249) printf("\n --> %lf %lf", uPrep[j], DummyY[j]>0?(1/(DummyY[j]))* (uPrep[j]) *(CPrep  - DummyC[i*3*MD->NumEle + j]):0.0);
//		if(j==0) printf("\n%lf %lf %lf %lf", DummyY[i], uPrep[j], CPrep, DummyC[i*3*MD->NumEle + j]);
		if(isnan(DummyDC[i*3*MD->NumEle + j]) == 1){ printf("debug 12 %d %d %lf %lf %lf %lf\n", i, j, DummyDC[i*3*MD->NumEle + j], MD->DummyY[j], uInfilUsat[j], uInfilSat[j]);}
		if(uInfilUsat[j]< 0 || uInfilSat[j] < 0)
		DummyDC[i*3*MD->NumEle + j] += MD->DummyY[j]>0.0?((-uInfilUsat[j]-uInfilSat[j])/(-uInfilUsat[j]-uInfilSat[j]+MD->DummyY[j])) * (CInfil - DummyC[i*3*MD->NumEle + j]):0.0;//TODO
		AGE DummyAgeDC[i*3*MD->NumEle + j] += MD->DummyY[j]>0.0?((-uInfilUsat[j]-uInfilSat[j])/(-uInfilUsat[j]-uInfilSat[j]+MD->DummyY[j])) * (AgeCInfil - DummyAgeC[i*3*MD->NumEle + j]):0.0;//TODO
if(fabs(DummyDC[i*3*MD->NumEle + j]) > 5*24*60){ printf("E$18 %d",j); getchar;}

		if(uInfilUsat[j] > 0)
		DummyDC[i*3*MD->NumEle + 1*MD->NumEle + j] += MD->DummyY[1*MD->NumEle + j]>0?(uInfilUsat[j]/(uInfilUsat[j]+MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity)) * (CInfil-DummyC[i*3*MD->NumEle + 1*MD->NumEle + j]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 1*MD->NumEle + j]) > 5*24*60){ printf("E$19 %d",j); getchar;}
		AGE DummyAgeDC[i*3*MD->NumEle + 1*MD->NumEle + j] += MD->DummyY[1*MD->NumEle + j]>0?(uInfilUsat[j]/(uInfilUsat[j]+MD->DummyY[1*MD->NumEle + j]*MD->Ele[j].Porosity)) * (AgeCInfil-DummyAgeC[i*3*MD->NumEle + 1*MD->NumEle + j]):0.0;
	
		if(uInfilSat[j] > 0)
		DummyDC[i*3*MD->NumEle + 2*MD->NumEle + j] += MD->DummyY[2*MD->NumEle + j]>0?(uInfilSat[j]/(uInfilSat[j]+MD->DummyY[2*MD->NumEle + j]*MD->Ele[j].Porosity)) * (CInfil-DummyC[i*3*MD->NumEle + 2*MD->NumEle + j]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + j]) > 5*24*60){ printf("E$20 %d",j); getchar;}
		AGE DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle + j] += MD->DummyY[2*MD->NumEle + j]>0?(uInfilSat[j]/(uInfilSat[j]+MD->DummyY[2*MD->NumEle + j]*MD->Ele[j].Porosity)) * (AgeCInfil-DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + j]):0.0;
		//?? DummyDC[i*3*MD->NumEle + j] += (1/(DummyY[j]+ABS_TOL))*(-5*uPrep[i])*(0      - DummyC[i*3*MD->NumEle + j]);
		if(isnan(DummyDC[i*3*MD->NumEle + j]) == 1){ printf("debug 13 %d %d %lf %lf %lf %lf\n", i, j, DummyDC[i*3*MD->NumEle + j], MD->DummyY[j], uInfilUsat[i], uInfilSat[i]);}

		//printf("%lf\t", uETSurf[j]);
		DummyDC[i*3*MD->NumEle + j] += (-uETSurf[j]+MD->DummyY[j])>0?(-uETSurf[j]/(-uETSurf[j]+MD->DummyY[j])) * (0 - DummyC[i*3*MD->NumEle + j]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + j]) > 5*24*60){ printf("E$21 %d",j); getchar;}
		//??DummyDC[i*3*MD->NumEle + 1*MD->NumEle + j] += (-uETUsat[j]+MD->DummyY[MD->NumEle + j])>0?(-uETUsat[j]/(-uETUsat[j]+MD->DummyY[MD->NumEle + j]*MD->Ele[j].Porosity)) * (0 - DummyC[i*3*MD->NumEle + 1*MD->NumEle + j]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 1*MD->NumEle + j]) > 5*24*60){ printf("E$22 %d",j); getchar;}
		AGE DummyAgeDC[i*3*MD->NumEle + j] += (-uETSurf[j]+MD->DummyY[j])>0?(-uETSurf[j]/(-uETSurf[j]+MD->DummyY[j])) * (0 - DummyAgeC[i*3*MD->NumEle + j]):0.0;

		//?? TODO add if statement DummyDC[i*3*MD->NumEle + 2*MD->NumEle + j] += (-uETSat[j]+MD->DummyY[2*MD->NumEle + j])>0?(-uETSat[j]/(-uETSat[j]+MD->DummyY[2*MD->NumEle + j]*MD->Ele[j].Porosity)) * (0 - DummyC[i*3*MD->NumEle + 2*MD->NumEle + j]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + j]) > 5*24*60){ printf("E$23 %d",j); getchar;}
		//?? TODO add if statement DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle + j] += (-uETSat[j]+MD->DummyY[2*MD->NumEle + j])>0?(-uETSat[j]/(-uETSat[j]+MD->DummyY[2*MD->NumEle + j]*MD->Ele[j].Porosity)) * (0 - DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + j]):0.0;

	}
    }




/**********************************************************************************/

// * RIVER-BED COMPONENT ONLY : START * /




// * RIVER-BED 2 RIVER-BED: START * /
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumRiv; j++){
			Ab =(MD->Riv[j].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,MD->DummyY[j + MD->NumRiv + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,3));

			// * RIVER & RIVER + DOWN : START * /
			if(qRiv[j][10]>0){
				Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
				AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
			}				
			else{
				if(MD->Riv[j].down>0){
					Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down -1];
					AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down -1];
				}
				else{
					//BC TODO
				}
			}
		
			DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+MD->NumRiv+j]>0?((-qRiv[j][10])/(-qRiv[j][10]+Ab*MD->DummyY[3*MD->NumEle+MD->NumRiv+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j ]) > 5*24*60){ printf("E$5 %d",j); getchar;}
			AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+ MD->NumRiv +j]>0?((-qRiv[j][10])/(-qRiv[j][10]+Ab*MD->DummyY[3*MD->NumEle+ MD->NumRiv +j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;

			if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]) == 1){ printf("debug 5 %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]);}
			
			if(MD->Riv[j].down>0){
				Ab_down =(MD->Riv[MD->Riv[j].down-1].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].interpOrd,MD->DummyY[MD->Riv[j].down-1 + 3*MD->NumEle + MD->NumRiv],MD->Riv_Shape[MD->Riv[MD->Riv[j].down-1].shape - 1].coeff,3));

				DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down-1] += MD->DummyY[3*MD->NumEle+ MD->NumRiv +MD->Riv[j].down-1]>0?((qRiv[j][10])/(qRiv[j][10]+Ab_down*MD->DummyY[3*MD->NumEle+ MD->NumRiv +MD->Riv[j].down-1]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down-1]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down-1 ]) > 5*24*60){ printf("E$6 %d",j); getchar;}
				AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down-1] += MD->DummyY[3*MD->NumEle+ MD->NumRiv +MD->Riv[j].down-1]>0?((qRiv[j][10])/(qRiv[j][10]+Ab_down*MD->DummyY[3*MD->NumEle+ MD->NumRiv+MD->Riv[j].down-1]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down-1]):0.0;

				if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down-1]) == 1){ printf("debug 6 %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + MD->Riv[j].down-1]);}
			}
		}
	}
// * RIVER-BED 2 RIVER-BED: END   * /







// * RIVER-BED 2 RIVER: START * /
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumRiv; j++){
			Ab =(MD->Riv[j].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,MD->DummyY[j + MD->NumRiv + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,3));

			if(qRiv[j][6]>0){
				Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
				AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
			}				
			else{
				Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
				AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
			}

			DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+ j]>0?((-qRiv[j][6])/(-qRiv[j][6]+Ab*MD->DummyY[3*MD->NumEle+ j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;

			DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+ MD->NumRiv + j]>0?((qRiv[j][6])/(qRiv[j][6]+Ab*MD->DummyY[3*MD->NumEle+ MD->NumRiv + j]*MD->Ele[j+MD->NumEle].Porosity))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;

			AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+ j]>0?((-qRiv[j][6])/(-qRiv[j][6]+Ab*MD->DummyY[3*MD->NumEle+ j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]):0.0;

			AGE DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+ MD->NumRiv + j]>0?((qRiv[j][6])/(qRiv[j][6]+Ab*MD->DummyY[3*MD->NumEle+ MD->NumRiv + j]*MD->Ele[j+MD->NumEle].Porosity))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;
		}
	}
// * RIVER-BED 2 RIVER: END   * /






// * RIVER-BED 2 SUB-SURFACE: START * /
	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<MD->NumRiv; j++){
			Ab =(MD->Riv[j].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,MD->DummyY[j + MD->NumRiv + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[j].shape - 1].coeff,3));

			if(qRiv[j][7] > 0){
                                Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
                                AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
			}				
                        else{
                                Cs = DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle - 1];
                                AgeCs = DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle - 1];
			}
                        DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+MD->NumRiv+j]>0?((-qRiv[j][7])/(-qRiv[j][7]+Ab*MD->DummyY[3*MD->NumEle+MD->NumRiv+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j ]) > 5*24*60){ printf("E$31 %d",j); getchar;}
                        DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]>0?((qRiv[j][7])/(qRiv[j][7]+MD->Ele[MD->Riv[j].LeftEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]*MD->Ele[MD->Riv[j].LeftEle-1].Porosity))*(Cs-DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1 ]) > 5*24*60){ printf("E$32 %d",j); getchar;}
			if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]) == 1){ printf("debug 7bb %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]);}                         
			if(isnan(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]) == 1){ printf("debug 8bb %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]);}
		AGE{
			DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+MD->NumRiv+j]>0?((-qRiv[j][7])/(-qRiv[j][7]+Ab*MD->DummyY[3*MD->NumEle+MD->NumRiv+j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;
                        DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]>0?((qRiv[j][7])/(qRiv[j][7]+MD->Ele[MD->Riv[j].LeftEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].LeftEle-1]*MD->Ele[MD->Riv[j].LeftEle-1].Porosity))*(AgeCs-DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].LeftEle-1]):0.0;
		}

                        if(qRiv[j][8] > 0){
                                Cs = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
                                AgeCs = DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j];
			}				
                        else{
                                Cs = DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle - 1];
                                AgeCs = DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle - 1];
			}				
                        DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+MD->NumRiv+j]>0?((-qRiv[j][8])/(-qRiv[j][8]+Ab*MD->DummyY[3*MD->NumEle+MD->NumRiv+j]))*(Cs-DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;
if(fabs(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]) > 5*24*60){ printf("E$33 %d",j); getchar;}
                        DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]>0?((qRiv[j][8])/(qRiv[j][8]+MD->Ele[MD->Riv[j].RightEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]*MD->Ele[MD->Riv[j].RightEle-1].Porosity))*(Cs-DummyC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]):0.0;
if(fabs(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]) > 5*24*60){ printf("E$34 %d",j); getchar;}
			if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]) == 1){ printf("debug 7cc %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]);}
			if(isnan(DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]) == 1){ printf("debug 8cc %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]);}
		AGE{
			DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j] += MD->DummyY[3*MD->NumEle+MD->NumRiv+j]>0?((-qRiv[j][8])/(-qRiv[j][8]+Ab*MD->DummyY[3*MD->NumEle+MD->NumRiv+j]))*(AgeCs-DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + MD->NumRiv + j]):0.0;
                        DummyAgeDC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1] += MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]>0?((qRiv[j][8])/(qRiv[j][8]+MD->Ele[MD->Riv[j].RightEle-1].area*MD->DummyY[2*MD->NumEle + MD->Riv[j].RightEle-1]*MD->Ele[MD->Riv[j].RightEle-1].Porosity))*(AgeCs-DummyAgeC[i*3*MD->NumEle + 2*MD->NumEle + MD->Riv[j].RightEle-1]):0.0;
		}
	     }
      	}
// * RIVER-BED 2 SUB-SURFACE: END   * /






// * RIVER-BED COMPONENT ONLY : END * /



               // for (i =0;i<10;i++){
               // printf("%lf\n ",NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv) + 0*3*MD->NumEle + 2*MD->NumEle + i));
               // }

//printf("%lf\t%lf\n",Y[(3*MD->NumEle+2*MD->NumRiv) + 0*3*MD->NumEle + 2*MD->NumEle + cid],DummyC[0*3*MD->NumEle + 2*MD->NumEle + cid]);

// * Assignment DummyDC -> DY : START * /
    for(i=0; i<MD->NumSolute; i++){
	for(j=0; j<3*MD->NumEle; j++){
		if(isnan(DummyDC[i*3*MD->NumEle + j]) == 1 || isinf(DummyDC[i*3*MD->NumEle + j]) == 1){ printf("ELE FLAG %d %d %lf\n", i, j, DummyDC[i*3*MD->NumEle + j]); getchar(); getchar();}
		//TODO if(DummyDC[i*3*MD->NumEle + j] < 0.0 && Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] <= ABS_TOL){ DummyDC[i*3*MD->NumEle + j] = 0.0; printf("getcharA1"); getchar();}
		if(DummyDC[i*3*MD->NumEle + j] < 0.0 && Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] <= 0){ DummyDC[i*3*MD->NumEle + j] = 0.0; printf("%d \t %lf \t %lf getcharA1",j,DummyDC[i*3*MD->NumEle + j], Y[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j]); getchar();}
		//??if(Y[j] <= ABS_TOL) DummyDC[i*3*MD->NumEle + j] = 0.0;
		if(Y[j] <= 0) DummyDC[i*3*MD->NumEle + j] = 0.0;
		DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] = (DummyDC[i*3*MD->NumEle + j]);///(60.0*24.0));
              //  if (j == 2*MD->NumEle + cid){
              //  printf("%lf",DY[(3*MD->NumEle+2*MD->NumRiv) + 0*3*MD->NumEle + 2*MD->NumEle + cid]);
              //  getchar();
              //  }
		//DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] = DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j]*MD->DummyY[j]*MD->Ele[j%MD->NumEle].area;
		//DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] = ((DummyDC[i*3*MD->NumEle + j]/(60.0*24.0))*MD->DummyY[j]+(DY[j]/(60.0*24.0))*DummyC[i*3*MD->NumEle + j])*MD->Ele[j%MD->NumEle].area;
//DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j] = 0.0;//??

		//if(j>2*MD->NumEle) printf(" %lf~%lf ",DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j],DY[3*MD->NumEle + j]);

		//AGE if(DummyAgeDC[i*3*MD->NumEle + j] < 0.0 && Y[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j] <= ABS_TOL) DummyAgeDC[i*3*MD->NumEle + j] = 0.0;
		AGE if(DummyAgeDC[i*3*MD->NumEle + j] < 0.0 && Y[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j] <= 0) DummyAgeDC[i*3*MD->NumEle + j] = 0.0;
		//AGE if(Y[j] <= ABS_TOL) DummyAgeDC[i*3*MD->NumEle + j] = 0.0;
		AGE if(Y[j] <= 0) DummyAgeDC[i*3*MD->NumEle + j] = 0.0;
		AGE DY[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j] = (DummyC[i*3*MD->NumEle + j]/UNIT_C + DummyAgeDC[i*3*MD->NumEle + j]);///(60.0*24.0);
		//AGE DY[(3*MD->NumEle+MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j] = DummyAgeDC[i*3*MD->NumEle + j]/(60.0*24.0);
		//if(j==2844){printf(" %.10lf %.10lf %.10lf \n", Y[j], DY[j], DY[(3*MD->NumEle+MD->NumRiv) + i*3*MD->NumEle + j]);}
		
		//DY[(3*MD->NumEle+MD->NumRiv) + i*3*MD->NumEle + j] = 0;
	}//getchar();
    }

	for(i=0; i<MD->NumSolute; i++){
		for(j=0; j<2*MD->NumRiv; j++){
		
		if(isnan(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1 || isinf(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]) == 1){ printf("RIV FLAG %d %d %lf\n", i, j, DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]); printf("getcharR1"); getchar(); getchar();}

		// TODO if(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] < 0.0 && Y[(3*MD->NumEle+2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] <= ABS_TOL){ DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] = 0.0;
		if(DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] < 0.0 && Y[(3*MD->NumEle+2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] <= 0){ DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] = 0.0;
//??printf("getcharR2"); getchar(); 
}

		//if(Y[3*MD->NumEle+j] <= ABS_TOL){ DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] = 0.0; //printf("getcharR3"); getchar();
		if(Y[3*MD->NumEle+j] <= 0){ DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] = 0.0;
}

		DY[(3*MD->NumEle+2*MD->NumRiv) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = DummyDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];///(60.0*24.0);


		AGE{ if(DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] < 0.0 && Y[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] <= 0){ DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] = 0.0;} }

		AGE{ if(Y[3*MD->NumEle+j] <= 0){ DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j] = 0.0;} }

		AGE DY[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] = DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j]/UNIT_C + DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j];
		//#AGE { if(j==0*MD->NumRiv) printf("%lf %lf %lf %lf %lf", DummyC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j], DummyAgeDC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j], DY[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j], DummyAgeC[MD->NumSolute*3*MD->NumEle + i*2*MD->NumRiv + j], MD->RivAgeConc[j][i]); }

		}
	}

/* Assignment DummyDC -> DY : END */

//Yi Update new information to state variable 

      if (((int)t % 60)==0){
        fprintf(omf[0],"%f\t", t);
        fprintf(omf[1],"%f\t", t);    
        fprintf(omf[2],"%f\t", t);
        fprintf(omf[3],"%f\t", t);    
      }
  somsobT = 0.0;
  somsolT = 0.0;
for (i=0; i<MD->NumSolute; i++){
	for(j=0; j<MD->NumEle; j++){
            AquiferDepth=(MD->Ele[j].zmax-MD->Ele[j].zmin); //Yi
            Deficit = AquiferDepth - MD->gwtb_old[j];//Yi
            MD->EleConc[0*MD->NumEle + j][i] = DummyC[i*3*MD->NumEle + 0*MD->NumEle + j] + DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle +0*MD->NumEle + j];
            MD->EleConc[1*MD->NumEle + j][i] = DummyC[i*3*MD->NumEle + 1*MD->NumEle + j] + DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle +1*MD->NumEle + j];
            MD->EleConc[2*MD->NumEle + j][i] = DummyC[i*3*MD->NumEle + 2*MD->NumEle + j] + DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle +2*MD->NumEle + j];

            DummyDS[0*MD->NumEle + j] = (MD->DummyY[1*MD->NumEle + j]>0)? -(MD->Tau/UNIT_C)*(DummyS[0*MD->NumEle + j]-MD->Kd*DummyC[i*3*MD->NumEle + 1*MD->NumEle + j]):0.0;
            DummyDS[1*MD->NumEle + j] = (MD->DummyY[2*MD->NumEle + j]>0)? -(MD->Tau/UNIT_C)*(DummyS[1*MD->NumEle + j]-MD->Kd*DummyC[i*3*MD->NumEle + 2*MD->NumEle + j]):0.0;

	    if (Deficit >  MD->unsatY_old[j]){
               if (MD->gwtb_old[j] > 0){
                  xmas = fabs(MD->Ele[j].S0/MD->somAlpha*exp(-MD->somAlpha*AquiferDepth)*(exp(MD->somAlpha*MD->DummyY[2*MD->NumEle + j]) - exp(MD->somAlpha*MD->gwtb_old[j])))*MD->Bd*MD->Ele[j].area;
                  DummyDSmas[0*MD->NumEle + j] =  DummyDS[0*MD->NumEle + j]*MD->Bd*MD->Ele[j].area*(AquiferDepth - MD->gwtb_old[j]);
                  DummyDSmas[0*MD->NumEle + j] =  (MD->DummyY[2*MD->NumEle + j] > MD->gwtb_old[j])? DummyDSmas[0*MD->NumEle + j] - xmas:DummyDSmas[0*MD->NumEle + j] + xmas; 

                  DummyDSmas[1*MD->NumEle + j] =  DummyDS[1*MD->NumEle + j]*MD->Bd*MD->Ele[j].area*(MD->gwtb_old[j]);
                  DummyDSmas[1*MD->NumEle + j] =  (MD->DummyY[2*MD->NumEle + j] > MD->gwtb_old[j])? DummyDSmas[1*MD->NumEle + j] + xmas:DummyDSmas[1*MD->NumEle + j] - xmas;
               }else{
                  DummyDSmas[0*MD->NumEle + j] =  DummyDS[0*MD->NumEle + j]*MD->Bd*MD->Ele[j].area*AquiferDepth;
                  DummyDSmas[1*MD->NumEle + j] =  0.0;
               }
            } else{
               xmas = fabs(MD->Ele[j].S0/MD->somAlpha*(exp(-MD->somAlpha*MD->DummyY[1*MD->NumEle + j]) - exp(-MD->somAlpha*MD->unsatY_old[j])))*MD->Bd*MD->Ele[j].area;
               DummyDSmas[0*MD->NumEle + j] = MD->Bd*MD->Ele[j].area*(MD->unsatY_old[j])*DummyDS[0*MD->NumEle + j];
               DummyDSmas[0*MD->NumEle + j] = (MD->DummyY[1*MD->NumEle + j] > MD->unsatY_old[j])? DummyDSmas[0*MD->NumEle + j] + xmas:DummyDSmas[0*MD->NumEle + j] - xmas;

               DummyDSmas[1*MD->NumEle + j] = MD->Bd*MD->Ele[j].area*(AquiferDepth - MD->unsatY_old[j])*DummyDS[1*MD->NumEle + j];
               DummyDSmas[1*MD->NumEle + j] = (MD->DummyY[1*MD->NumEle + j] > MD->unsatY_old[j])? DummyDSmas[1*MD->NumEle + j] - xmas:DummyDSmas[1*MD->NumEle + j] + xmas;
            }
            MD->EleS[0*MD->NumEle + j] = DummySmas[0*MD->NumEle + j] + DummyDSmas[0*MD->NumEle + j];
            MD->EleS[1*MD->NumEle + j] = DummySmas[1*MD->NumEle + j] + DummyDSmas[1*MD->NumEle + j];
            
            Deficit = AquiferDepth - MD->DummyY[2*MD->NumEle + j];//Yi
	    if (Deficit >  MD->DummyY[1*MD->NumEle + j]){
            MD->EleSbar[0*MD->NumEle + j] = MD->DummyY[2*MD->NumEle + j]>0? MD->EleS[0*MD->NumEle + j]/(MD->Bd*MD->Ele[j].area*Deficit):MD->EleS[0*MD->NumEle + j]/(MD->Bd*MD->Ele[j].area*AquiferDepth);
            MD->EleSbar[1*MD->NumEle + j] = MD->DummyY[2*MD->NumEle + j]>0? MD->EleS[1*MD->NumEle + j]/(MD->Bd*MD->Ele[j].area*(MD->DummyY[2*MD->NumEle + j])):0.0;
            }else{
            MD->EleSbar[0*MD->NumEle + j] = MD->EleS[0*MD->NumEle + j]/(MD->Bd*MD->Ele[j].area*MD->DummyY[1*MD->NumEle + j]);
            MD->EleSbar[1*MD->NumEle + j] = MD->EleS[1*MD->NumEle + j]/(MD->Bd*MD->Ele[j].area*(AquiferDepth - MD->DummyY[1*MD->NumEle + j]));
            }
            somsobT = somsobT + MD->EleS[0*MD->NumEle + j]/1000000;
            somsobT = somsobT + MD->EleS[1*MD->NumEle + j]/1000000;
            somsolT = somsolT + MD->EleConc[1*MD->NumEle + j][i]*MD->Ele[j].area*MD->Ele[j].Porosity*MD->DummyY[1*MD->NumEle + j]/1000000;
            somsolT = somsolT + MD->EleConc[2*MD->NumEle + j][i]*MD->Ele[j].area*MD->Ele[j].Porosity*MD->DummyY[2*MD->NumEle + j]/1000000;
	    // Calculate the SOM profile for new time step
	    MD->Ele[j].S0 =(MD->EleS[0*MD->NumEle + j]+MD->EleS[1*MD->NumEle + j])*(MD->somAlpha)/(MD->Bd*MD->Ele[j].area*(1 - exp(-MD->somAlpha*(AquiferDepth))));
        if (((int)t % 60)==0){
            fprintf(omf[0],"%lf\t", MD->EleSbar[0*MD->NumEle + j]);
            fprintf(omf[1],"%lf\t", MD->EleSbar[1*MD->NumEle + j]);
            fprintf(omf[2],"%lf\t", MD->EleS[0*MD->NumEle + j]/1000000);
            fprintf(omf[3],"%lf\t", MD->EleS[1*MD->NumEle + j]/1000000);
           }
       }
}
        

        if (((int)t % 60)==0){
           fprintf(omf[0],"\n");
           fprintf(omf[1],"\n");
           fprintf(omf[2],"\n");
           fprintf(omf[3],"\n");
          }
    fprintf(omf[4],"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",MD->EleConc[1*MD->NumEle + cid][0],MD->EleConc[2*MD->NumEle + cid][0],MD->EleSbar[0*MD->NumEle + cid],MD->EleSbar[1*MD->NumEle + cid],MD->EleS[0*MD->NumEle + cid]/1000000,MD->EleS[1*MD->NumEle + cid]/1000000,somsolT,somsobT);
   fprintf(omf[4],"\n");
//Yi



//printf("\n#");

// * CALCULATION OF THE NEW SOLUTE CONCENTRATION: START */
for (i=0; i<MD->NumSolute; i++){
	for(j=0; j<3*MD->NumEle; j++){
//		NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv) + i*3*MD->NumEle + j) = DummyC[i*3*MD->NumEle + j] + DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j];
//		MD->EleConc[j][i]=DummyC[i*3*MD->NumEle + j] + DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j];
                NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv) + i*3*MD->NumEle + j) = MD->EleConc[j][i];
		AGE{
		NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j) = DummyAgeC[i*3*MD->NumEle + j] + DY[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j];
		MD->EleAgeConc[j][i]=DummyAgeC[i*3*MD->NumEle + j] + DY[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j];
		}


		/*if(j<MD->NumEle && uPrep[j]>0){
				NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv) + i*3*MD->NumEle + j) = CPrep;
				MD->EleConc[j][i]=CPrep;
		}*/
	}
	for(j=0; j<2*MD->NumRiv; j++){
		//NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle+j + i*2*MD->NumRiv) += DY[(3*MD->NumEle+2*MD->NumRiv) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j];
		//NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv) + MD->NumSolute*3*MD->NumEle+j + i*2*MD->NumRiv) += 0.1;
		NV_Ith_S(CV_Y,(3*MD->NumEle+2*MD->NumRiv) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j) = DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] + DY[(3*MD->NumEle+2*MD->NumRiv) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j];
		MD->RivConc[j][i]=DummyC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] + DY[(3*MD->NumEle+2*MD->NumRiv) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j];
		AGE{
		NV_Ith_S(CV_Y,(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j) = DummyAgeC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] + DY[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j];
		MD->RivAgeConc[j][i]=DummyAgeC[(MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j] + DY[(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j];
		}
		//printf("%lf\t",NV_Ith_S(CV_Y,(3*MD->NumEle+2*MD->NumRiv) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j));
	}

	// WORK IN PROGRESS
	for(j=0*MD->NumEle; j<3*MD->NumEle; j++){
		if(MD->DummyY[j] <= 0 && j > MD->NumEle){
			NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv) + i*3*MD->NumEle + j) = 0.0;
			MD->EleConc[j][i] = 0.0;
		}
		if(MD->DummyY[j] <= 0){
			AGE{
			NV_Ith_S(CV_Y, (3*MD->NumEle + 2*MD->NumRiv)*(MD->NumSolute+1) + i*3*MD->NumEle + j) = 0.0;
			MD->EleAgeConc[j][i] = 0.0;
			}

		}
		//if(j==2*MD->NumEle+2){ printf("| %lf %lf %lf || %lf %lf %lf || %lf %lf %lf # %lf %lf %lf\n", MD->DummyY[j-2*MD->NumEle], MD->EleConc[j-2*MD->NumEle][i], 1000000.0*DY[j-2*MD->NumEle], MD->DummyY[j-1*MD->NumEle], MD->EleConc[j-1*MD->NumEle][i], DY[j-1*MD->NumEle]==0?0.0:DY[j-1*MD->NumEle]>0?1.0:-1.0, MD->DummyY[j-0*MD->NumEle], MD->EleConc[j-0*MD->NumEle][i], 1000000.0*DY[j-0*MD->NumEle], DummyC[i*3*MD->NumEle + j], MD->Recharge[j-2*MD->NumEle], DY[(3*MD->NumEle+2*MD->NumRiv) + i*3*MD->NumEle + j]); }
	}
	for(j=1*MD->NumRiv; j<2*MD->NumRiv; j++){
		if(MD->DummyY[3*MD->NumEle+j] <= 0){
			NV_Ith_S(CV_Y,(3*MD->NumEle+2*MD->NumRiv) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j) = 0.0;
			MD->RivConc[j][i] = 0.0;

			AGE{
			NV_Ith_S(CV_Y,(3*MD->NumEle+2*MD->NumRiv)*(MD->NumSolute+1) + (MD->NumSolute*3*MD->NumEle) + i*2*MD->NumRiv + j) = 0.0;
			MD->RivAgeConc[j][i] = 0.0;
			}
		}
	}
}
// * CALCULATION OF THE NEW SOLUTE CONCENTRATION: END   */

//Yi Update gw table old at new time step

    for (i=0; i<MD->NumEle;i++){
//        printf("%lf\n", MD->gwtb_old[i] - MD->DummyY[i+2*MD->NumEle]);
        MD->gwtb_old[i] = (MD->DummyY[i+2*MD->NumEle]>0)? MD->DummyY[i+2*MD->NumEle] : 0.0 ; //Yi
        MD->unsatY_old[i] = MD->DummyY[i + 1*MD->NumEle] ;
        }
//Yi


    for(i=0; i<MD->NumEle; i++){
        free(uSub[i]);
        free(uSurf[i]);
        free(Disp[i]);
    }

    free(uSub); free(Disp);
    free(uSurf);

    for(i=0; i<MD->NumSolute; i++){
	free(Kd[i]); free(R[i]);
    }
    free(Kd); free(R);

    for(i=0; i<MD->NumRiv; i++){
        free(uRiv[i]);
        free(qRiv[i]);
    }
    free(uRiv); free(qRiv); free(uRivCount); free(DispRiv);

    free(uPrep); free(uInfilUsat); free(uInfilSat); free(uETSurf); free(uETUsat); free(uETSat);

    free(DummyC);  free(DummyAgeC);
    free(DummyDC); free(DummyAgeDC);
    free(DummyS); free(DummyDS); //Yi
    free(DummySmas); free(DummyDSmas); //Yi
    
    free(DY);


/* TRANSPORT */

	return 0;
	}

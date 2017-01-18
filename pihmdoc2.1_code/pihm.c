#define TRANSPORT
/*******************************************************************************
 * File        : pihm.c                                                        *
 * Version     : Nov, 2007 (2.0)                                               *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) All modifications in physical process representations  in this version   *
 *    are listed as header in f.c and is_sm_et.c.     			       *
 * b) All addition/modifications in variable and structure definition/declarat-*
 *    -ion are listed as header in read_alloc.c and initialize.c	       *
 * c) 3 new input files have been added for geology, landcover and calibration *
 *    data								       *
 * d) Ported to Sundials 2.1.0                                                 *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * PIHM is an integrated finite volume based hydrologic model. It simulates    * 
 * channel routing, overland flow, groundwater flow, macropore based infiltra- *
 * tion and stormflow, throughfall, evaporation from overlandflow-subsurface-  *
 * canopy, transpiration and  snowmelt by full coupling of processes.          * 
 * It uses semi-discrete finite volume approach to discretize PDEs into ODEs,  * 
 * and henceforth solving the global system of ODEs using CVODE. Global ODEs   *
 * are created in f.c. Any modifications in the process equations has to be    *
 * performed in f.c
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact                                   *
 *      --> Mukesh Kumar (muk139@psu.edu)                                      *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                  *
 * This code is free for research purpose only.                                *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *									       *
 * DEVELOPMENT RELATED REFERENCES:					       *
 * PIHM2.0:								       *
 *	a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *	Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *	b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *	Mesoscale Watershed", Advances in Water Resources (submitted)          *
 * PIHM1.0:								       *
 *	a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *	simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *	b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *	for multiprocess watershed simulation". Water Resources Research       *
 *-----------------------------------------------------------------------------*
 * LICENSE: 
 *******************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* SUNDIAL Header Files */
#include "sundials_types.h"   /* realtype, integertype, booleantype defination */
#include "cvode.h"           /* CVODE header file                             */
#include "cvode_spgmr.h"         /* CVSPGMR linear header file                    */
#include "sundials_smalldense.h"      /* use generic DENSE linear solver for "small"   */
#include "nvector_serial.h"  /* contains the definition of type N_Vector      */
#include "sundials_math.h"    /* contains UnitRoundoff, RSqrt, SQR functions   */
#include "cvode_dense.h"         /* CVDENSE header file                           */
#include "sundials_dense.h"           /* generic dense solver header file              */
#include "pihm.h"            /* Data Model and Variable Declarations     */
#include "f_sol.h"
#define UNIT_C 1440	     /* Unit Conversions */	

/* Function Declarations */
void initialize(char *, Model_Data, Control_Data *, N_Vector);
void is_sm_et(realtype, realtype, Model_Data, N_Vector);	
/* Function to calculate right hand side of ODE systems */
int f(realtype, N_Vector, N_Vector, void *);
void read_alloc(char *, Model_Data, Control_Data *);	/* Variable definition */
void update(realtype, Model_Data);	 
void PrintData(FILE **,Control_Data *, Model_Data, N_Vector, realtype);

/* Main Function */
int main(int argc, char *argv[])
	{  
	char tmpLName[20],tmpFName[20];	/* rivFlux File names */
  	Model_Data mData;               /* Model Data                */
  	Control_Data cData;             /* Solver Control Data       */
  	N_Vector CV_Y;                  /* State Variables Vector    */
  	void *cvode_mem;                /* Model Data Pointer        */
  	int flag;                       /* flag to test return value */
  	FILE *Ofile[33];           	/* Output file     */
	char *ofn[33];
	FILE *iproj;			/* Project File */
  	int N;                          /* Problem size              */
  	int i,j,k;                      /* loop index                */
  	realtype t;    			/* simulation time           */
  	realtype NextPtr, StepSize;     /* stress period & step size */
  	clock_t start, end_r, end_s;    /* system clock at points    */
  	realtype cputime_r, cputime_s;  /* for duration in realtype  */
	char *filename;
//Yi
        FILE *OMfile[6];
        char *omfn[6];
        FILE *masbal;  //Yi
//        realtype *gwtb_old;
	/* Project Input Name */
	if(argc!=2)
		{
		iproj=fopen("projectName.txt","r");
		if(iproj==NULL)
			{
			printf("\t\nUsage ./pihm project_name");
			printf("\t\n         OR              ");
			printf("\t\nUsage ./pihm, and have a file in the current directory named projectName.txt with the project name in it");
			exit(0);
			}
		else
			{
			filename = (char *)malloc(15*sizeof(char));
			fscanf(iproj,"%s",filename);
			}
		}
	else
		{
  		/* get user specified file name in command line */
    		filename = (char *)malloc(strlen(argv[1])*sizeof(char));
		strcpy(filename,argv[1]);
		}
	/* Open Output Files */
	ofn[0] = (char *)malloc((strlen(filename)+3)*sizeof(char));
	strcpy(ofn[0], filename);
	Ofile[0]=fopen(strcat(ofn[0], ".GW"),"w");
	ofn[1] = (char *)malloc((strlen(filename)+5)*sizeof(char));
	strcpy(ofn[1], filename);
	Ofile[1]=fopen(strcat(ofn[1], ".surf"),"w");
	ofn[2] = (char *)malloc((strlen(filename)+4)*sizeof(char));
	strcpy(ofn[2], filename);
	Ofile[2]=fopen(strcat(ofn[2], ".et0"),"w");
	ofn[3] = (char *)malloc((strlen(filename)+4)*sizeof(char));
	strcpy(ofn[3], filename);
	Ofile[3]=fopen(strcat(ofn[3], ".et1"),"w");
	ofn[4] = (char *)malloc((strlen(filename)+4)*sizeof(char));
	strcpy(ofn[4], filename);
	Ofile[4]=fopen(strcat(ofn[4], ".et2"),"w");
	ofn[5] = (char *)malloc((strlen(filename)+3)*sizeof(char));
	strcpy(ofn[5], filename);
	Ofile[5]=fopen(strcat(ofn[5], ".is"),"w");
	ofn[6] = (char *)malloc((strlen(filename)+5)*sizeof(char));
	strcpy(ofn[6], filename);
	Ofile[6]=fopen(strcat(ofn[6], ".snow"),"w");
	for(i=0;i<11;i++)
		{
		sprintf(tmpLName,".rivFlx%d",i);
		strcpy(tmpFName,filename);
		strcat(tmpFName,tmpLName);
		Ofile[7+i]=fopen(tmpFName,"w");
		}	
	ofn[18] = (char *)malloc((strlen(filename)+6)*sizeof(char));
	strcpy(ofn[18], filename);
	Ofile[18]=fopen(strcat(ofn[18], ".stage"),"w");
	ofn[19] = (char *)malloc((strlen(filename)+6)*sizeof(char));
	strcpy(ofn[19], filename);
	Ofile[19]=fopen(strcat(ofn[19], ".unsat"),"w");
	ofn[20] = (char *)malloc((strlen(filename)+5)*sizeof(char));
	strcpy(ofn[20], filename);
	Ofile[20]=fopen(strcat(ofn[20], ".Rech"),"w");
	ofn[21] = (char *)malloc((strlen(filename)+5)*sizeof(char));//??BHATT +3
        strcpy(ofn[21], filename);
        Ofile[21]=fopen(strcat(ofn[21], ".rbed"),"w");
	ofn[22] = (char *)malloc((strlen(filename)+5)*sizeof(char));//??BHATT +3
        strcpy(ofn[22], filename);
        Ofile[22]=fopen(strcat(ofn[22], ".infil"),"w");
	
	ofn[23] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        strcpy(ofn[23], filename);
        Ofile[23]=fopen(strcat(ofn[23], ".GWsol"),"w");
	ofn[24] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        strcpy(ofn[24], filename);
        Ofile[24]=fopen(strcat(ofn[24], ".USsol"),"w");
	ofn[25] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        strcpy(ofn[25], filename);
        Ofile[25]=fopen(strcat(ofn[25], ".OLsol"),"w");

	ofn[26] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        strcpy(ofn[26], filename);
        Ofile[26]=fopen(strcat(ofn[26], ".RSsol"),"w");
	ofn[27] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        strcpy(ofn[27], filename);
        Ofile[27]=fopen(strcat(ofn[27], ".RBsol"),"w");


	//ofn[28] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        //strcpy(ofn[28], filename);
        //Ofile[28]=fopen(strcat(ofn[28], ".GWage"),"w");
	//ofn[29] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        //strcpy(ofn[29], filename);
        //Ofile[29]=fopen(strcat(ofn[29], ".USage"),"w");
	//ofn[30] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        //strcpy(ofn[30], filename);
        //Ofile[30]=fopen(strcat(ofn[30], ".OLage"),"w");

	//ofn[31] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        //strcpy(ofn[31], filename);
        //Ofile[31]=fopen(strcat(ofn[31], ".RSage"),"w");
	//ofn[32] = (char *)malloc((strlen(filename)+6)*sizeof(char));//??BHATT +3
        //strcpy(ofn[32], filename);
        //Ofile[32]=fopen(strcat(ofn[32], ".RBage"),"w");

//Yi 
        omfn[0] = (char *)malloc((strlen(filename)+6)*sizeof(char));
        strcpy(omfn[0],filename);
        OMfile[0] = fopen(strcat(omfn[0], ".USsom"),"w");
        omfn[1] = (char *)malloc((strlen(filename)+6)*sizeof(char));
        strcpy(omfn[1],filename);
        OMfile[1] = fopen(strcat(omfn[1], ".GWsom"),"w");
        omfn[2] = (char *)malloc((strlen(filename)+6)*sizeof(char));
        strcpy(omfn[2],filename);
        OMfile[2] = fopen(strcat(omfn[2], ".USmas"),"w");
        omfn[3] = (char *)malloc((strlen(filename)+6)*sizeof(char));
        strcpy(omfn[3],filename);
        OMfile[3] = fopen(strcat(omfn[3], ".GWmas"),"w");
        omfn[4] = (char *)malloc((strlen(filename)+6)*sizeof(char));
        strcpy(omfn[4],filename);
        OMfile[4] = fopen(strcat(omfn[4], ".Einfo"),"w");
        omfn[5] = (char *)malloc((strlen(filename)+7)*sizeof(char));
        strcpy(omfn[5],filename);
        OMfile[5] = fopen(strcat(omfn[5], ".EinfoY"),"w");
         
        
  	/* allocate memory for model data structure */
  	mData = (Model_Data)malloc(sizeof *mData);
  
  	printf("\n ...  PIHM 2.0 is starting ... \n");
 
 	/* read in 9 input files with "filename" as prefix */
  	read_alloc(filename, mData, &cData); 

/*	if(mData->UnsatMode ==1)
		{    
  		} */
	if(mData->UnsatMode ==2)
		{    
  		/* problem size */
  		N = 3*mData->NumEle + 2*mData->NumRiv;
		TRANSPORT N = (3*mData->NumEle + 2*mData->NumRiv) * (1+mData->NumSolute*(1+mData->AGE));
		mData->DummyY=(realtype *)malloc((3*mData->NumEle+2*mData->NumRiv)*sizeof(realtype));
  		}
	printf("\nProblem Size N= %d\n", N);  
	
  	/* initial state variable depending on machine*/
  	CV_Y = N_VNew_Serial(N);
  
  	/* initialize mode data structure */
  	initialize(filename, mData, &cData, CV_Y);
  
  	printf("\nSolving ODE system ... \n");
  
  	/* allocate memory for solver */
  	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  	if(cvode_mem == NULL) {printf("CVodeMalloc failed. \n"); return(1);}
  
  	flag = CVodeSetFdata(cvode_mem, mData);  
  	flag = CVodeSetInitStep(cvode_mem,cData.InitStep);
  	flag = CVodeSetStabLimDet(cvode_mem,TRUE);  
  	flag = CVodeSetMaxStep(cvode_mem,cData.MaxStep); //?? COMMENT THIS FOR BETTER SIMULATION EFFICIENCY
  	flag = CVodeMalloc(cvode_mem, f, cData.StartTime, CV_Y, CV_SS, cData.reltol, &cData.abstol);  
  	flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
  	//flag = CVSpgmrSetGSType(cvode_mem, MODIFIED_GS);
  
  	/* set start time */
  	t = cData.StartTime;
  	start = clock();
	PrintData(Ofile,&cData,mData, CV_Y,t);
masbal = fopen("mas.chek","w"); //Yi
  	/* start solver in loops */
  	for(i=0; i<cData.NumSteps; i++)
  		{
	/*	if (cData.Verbose != 1)
    			{
      			printf("  Running: %-4.1f%% ... ", (100*(i+1)/((realtype) cData.NumSteps))); 
      			fflush(stdout);
    			} */
    		/* inner loops to next output points with ET step size control */
    		while(t < cData.Tout[i+1])
    			{
      			if (t + cData.ETStep >= cData.Tout[i+1])
      				{
        			NextPtr = cData.Tout[i+1];
      				}
      			else
      				{
        			NextPtr = t + cData.ETStep;
      				}
      			StepSize = NextPtr - t; 
      
      			/* calculate Interception Storage */
      			is_sm_et(t, StepSize, mData,CV_Y);
			printf("\n Tsteps = %f ",t);
      			flag = CVode(cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL); 
			update(t,mData);
    			}
//Yi for water mass balance
                for(k = 0; k < mData->NumEle; k++){
                        mData->tov = mData->tov + mData->Ele[k].area*NV_Ith_S(CV_Y,k + 0 * mData->NumEle)/mData->Tarea; //Yi
                        mData->tus = mData->tus + mData->Ele[k].Porosity*mData->Ele[k].area*NV_Ith_S(CV_Y, k + 1 * mData->NumEle)/mData->Tarea; //Yi
                        mData->tgw = mData->tgw + mData->Ele[k].Porosity*mData->Ele[k].area*NV_Ith_S(CV_Y, k + 2 * mData->NumEle)/mData->Tarea; //Yi
                        mData->tsnowc = mData->tsnowc + mData->Ele[k].area*mData->EleSnowCanopy[k]/mData->Tarea; //Yi
                        mData->tsnowg = mData->tsnowg + mData->Ele[k].area*mData->EleSnowGrnd[k]/mData->Tarea; //Yi
                        mData->tis = mData->tis + mData->Ele[k].area*mData->EleIS[k]/mData->Tarea; //Yi
                        mData->et0flx = mData->et0flx + mData->EleET[k][0]/1440*mData->Ele[k].area/mData->Tarea;
                        mData->et1flx = mData->et1flx + mData->EleET[k][1]/1440*mData->Ele[k].area/mData->Tarea;
                        mData->et2flx = mData->et2flx + mData->EleET[k][2]/1440*mData->Ele[k].area/mData->Tarea;
                        mData->pflx = mData->pflx + mData->ElePrep[k]/1440*mData->Ele[k].area/mData->Tarea;

                }
                for (k= 0; k < mData->NumRiv; k++ ){
                    mData->triv = mData->triv + mData->Riv[k].Length*CS_AreaOrPerem(mData->Riv_Shape[mData->Riv[k].shape - 1].interpOrd,NV_Ith_S(CV_Y, k + 3 * mData->NumEle), mData->Riv[k].coeff, 1)/mData->Tarea;   //Yi
                    mData->tbriv = mData->tbriv + NV_Ith_S(CV_Y, k + mData->NumRiv + 3 * mData->NumEle) * mData->Riv[k].Length*CS_AreaOrPerem(mData->Riv_Shape[mData->Riv[k].shape - 1].interpOrd, mData->Riv[k].depth, mData->Riv[k].coeff, 3)/mData->Tarea; //Yi
                }

                mData->rsflx = mData->rsflx + mData->FluxRiv[5][1]/1440/mData->Tarea;
                //fprintf(masbal, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",mData->inius,mData->inigw,mData->tus,mData->tgw, mData->et0flx, mData->et1flx, mData->et2flx, mData->rsflx,mData->pflx);//Yi
                //fprintf(masbal, "%lf %lf %lf\n",t,(mData->tov + mData->tgw + mData->tus + mData->tsnowc + mData->tsnowg + mData->tis) - (mData->iniov + mData->inius + mData->inigw), mData->pflx-(mData->et0flx + mData->et1flx + mData->et2flx + mData->rsflx));//Yi
                fprintf(masbal, "%lf %lf %lf %lf\n",t,(mData->tov + mData->tgw + mData->tus + mData->tsnowc + mData->tsnowg + mData->tis) - (mData->iniov + mData->inius + mData->inigw), mData->pflx-(mData->et0flx + mData->et1flx + mData->et2flx + mData->rsflx),(mData->triv+mData->tbriv)-(mData->iniriv+mData->inibriv));//Yi

                mData->tov = 0.0;
                mData->tus = 0.0;
                mData->tgw = 0.0;
                mData->tsnowc = 0.0;
                mData->tsnowg = 0.0;
                mData->tis = 0.0;
                mData->triv = 0.0;
                mData->tbriv = 0.0;
//Yi

		f_sol(OMfile, t, CV_Y, NULL, mData);
		PrintData(Ofile,&cData,mData, CV_Y,t);
  		}
fclose(masbal); //Yi
   	/* Free memory */
  	N_VDestroy_Serial(CV_Y);
  	/* Free integrator memory */
  	CVodeFree(&cvode_mem);
  	free(mData);
  	return 0;
	}


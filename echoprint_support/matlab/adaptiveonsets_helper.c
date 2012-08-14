/*
 * adaptiveonsets_helper.c
 * 
 * Matlab MEX implementation of core of adaptiveonsets.m
 *
 * Dan Ellis dpwe@ee.columbia.edu 2010-10-25
 */

#include "mex.h"
#include "string.h"
#include <math.h>

#ifndef MWSIZE_MIN
/* matlab R17 -ism? */
typedef unsigned int mwSize;
#endif

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int 	i,j,k;

    if (nrhs < 2){
	mexPrintf("[O,D,U,M] = adaptiveonsets_helper(E,T,DT)\n");
	mexPrintf("           E is a sgram-like matrix of energies.  Return in O a conformal\n");
	mexPrintf("           mask of where onsets are found.  D returns the decaying \n");
	mexPrintf("           threshold, and U the decay time. T defines target # steps\n");
	mexPrintf("           between successive hits.  M gives prominences.\n");
	mexPrintf("           DT defines deadtime, min gap between hits (dflt 32).\n");
	return;
    }

#define MAXROWS 65

    if (nlhs > 0){
	mxArray  *Omatrix, *Dmatrix, *Umatrix, *Mmatrix;
	double *pE, *pD, *pO, *pU, *pM;
	int rows, cols, i, j;
	int deadtime = 32;
	int ttarg;
	mwSize siz;
	double H[MAXROWS],taus[MAXROWS], N[MAXROWS], Y0[MAXROWS];
	int contact[MAXROWS], lcontact[MAXROWS], tsince[MAXROWS];
	//double overfact = 1.05;  /* threshold rel. to actual peak */
	double overfact = 1.02;  /* threshold rel. to actual peak */
	double bn[] = {0.1883, 0.4230, 0.3392}; /* preemph filter */
	int nbn = 3;
	double a1 = 0.98;

	rows = mxGetM(prhs[0]);
	if (rows > MAXROWS) {
	    mexPrintf("error: rows %d is larger than %d\n", rows, MAXROWS);
	    return;
	}

	cols = mxGetN(prhs[0]);
	pE = (double *)mxGetData(prhs[0]);
	
	ttarg = (int)rint(mxGetScalar(prhs[1]));

	if (nrhs > 2) {
	    deadtime = (int)rint(mxGetScalar(prhs[2]));
	}
	/* mexPrintf("nrhs = %d, deadtime = %d\n", nrhs, deadtime); */

 	/* mexPrintf("rows=%d cols=%d type=%s\n", rows, cols, mxGetClassName(prhs[0])); */

	Omatrix = mxCreateNumericMatrix(rows, cols, mxDOUBLE_CLASS, mxREAL);
	Dmatrix = mxCreateNumericMatrix(rows, cols, mxDOUBLE_CLASS, mxREAL);
	Umatrix = mxCreateNumericMatrix(rows, cols, mxDOUBLE_CLASS, mxREAL);
	Mmatrix = mxCreateNumericMatrix(rows, cols, mxDOUBLE_CLASS, mxREAL);

	pO = (double *)mxGetData(Omatrix);
	pD = (double *)mxGetData(Dmatrix);
	pU = (double *)mxGetData(Umatrix);
	pM = (double *)mxGetData(Mmatrix);

	for (j = 0; j < rows; ++j) {
	    taus[j] = 1.0;
	    H[j] = pE[j];
	    contact[j] = 0;
	    lcontact[j] = 0;
	    tsince[j] = 0;
	    Y0[j] = 0;
	}

	for (i = 0; i < cols; ++i) {
	    for (j = 0; j < rows; ++j) {
		double xn = 0;
		/* calculate the filter */
		/* FIR part */
		if (i >= 2*nbn) {
		    for (k = 0; k < nbn; ++k) {
			xn += bn[k]*(pE[j-rows*k] - pE[j-rows*(2*nbn-k)]);
		    }
		}
		/* IIR part */
		xn = xn + a1*Y0[j];

		/* defeat - make it like it used to be */
		/* xn = pE[j]; /* */

		/* apply threshold on output */
		contact[j] = (xn > H[j])? 1 : 0;
		if (contact[j] == 1 && lcontact[j] == 0) {
		    /* attach - record the threshold level unless we have one */
		    if(N[j] == 0) {
			N[j] = H[j];
		    }
		}
		if (contact[j] == 1) {
		    /* update with new threshold */
		    H[j] = xn * overfact;
		} else {
		    /* apply decays */
		    H[j] = H[j] * exp(-1.0/(double)taus[j]);
		}
		if (contact[j] == 0 && lcontact[j] == 1) {
		    /* detach */
		    pO[j] = 1;
		    pM[j] = Y0[j]/N[j];
		    if (pM[j] > 20.0) {
			pM[j] = 20.0;
		    }
		    /* apply deadtime */
		    for(k = 1; k < ((i > deadtime)?deadtime:i); ++k) {
			pO[j - k*rows] = 0;
			pM[j - k*rows] = 0;
		    }
		    tsince[j] = 0;
		}
		++tsince[j];
		if (tsince[j] > ttarg) {
		    taus[j] = taus[j] - 1;
		    if (taus[j] < 1) taus[j] = 1;
		} else {
		    taus[j] = taus[j] + 1;
		}
		if ( (contact[j] == 0) &&  (tsince[j] > deadtime)) {
		    /* forget the threshold where we recently hit */
		    N[j] = 0;
		}
		lcontact[j] = contact[j];
		//pD[j] = H[j];
		pD[j] = N[j];
		pU[j] = taus[j];
		/* remember the last filtered level */
		Y0[j] = xn;
	    }
	    pE += rows;
	    pO += rows;
	    pD += rows;
	    pU += rows;
	    pM += rows;
	}


	plhs[0] = Omatrix;
	if (nlhs > 1) {  plhs[1] = Dmatrix; }
	if (nlhs > 2) {  plhs[2] = Umatrix; }
	if (nlhs > 3) {  plhs[3] = Mmatrix; }
    }
}


/*
 * online_whiten_helper.c
 * 
 * Matlab MEX implementation of core of online_whiten.m
 *
 * Dan Ellis dpwe@ee.columbia.edu 2011-06-13
 */

#include "mex.h"
#include "string.h"
#include <math.h>

#ifndef MWSIZE_MIN
/* matlab R17 -ism? */
typedef unsigned int mwSize;
#endif

static float *R = NULL;
static float *Xo = NULL;
static float *ai = NULL;
static float *acc = NULL;
static int P = 0;
static float T = 0;


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int 	i,j,k;

    if (nrhs == 0){
	mexPrintf("[Y] = online_whiten_helper(X,P,T)\n");
	mexPrintf("           Stateful online whitening of X.\n");
	mexPrintf("           X is the (next) vector of input samples, P is LPC order, \n");
	mexPrintf("           T is time constant of autoco smoothing.\n");
	mexPrintf("           Pass in zero-length X to clear internal filter state\n");
	mexPrintf("           and accumulated autocorrelation\n");
	return;
    }

    if (nrhs > 0){
	mxArray  *Ymatrix;
	double *pX, *pY;
	int i, j, k, N;
	float alpha, E, ki;

	if (nrhs > 1) {
	    P = (int)(*(double *)mxGetData(prhs[1]));
	}
	if (nrhs > 2) {
	    T = (float)(*(double *)mxGetData(prhs[2]));
	}
	
	/* enforced defaults */
	if (P <= 0) { P = 40; }
	if (T <= 0) { T = 100; }

	alpha = 1.0/T;

	/* how many points passed in? */
	N = mxGetM(prhs[0]) * mxGetN(prhs[0]);
	pX = mxGetData(prhs[0]);

	/* special case to clear state */
	if (N == 0) {
	    /* mexPrintf("freeing stored tables\n"); */
	    if (R != NULL)  { free(R); R = NULL; }
	    if (Xo != NULL) { free(Xo); Xo = NULL; }
	    if (ai != NULL) { free(ai); ai = NULL; }
	    if (acc != NULL) { free(acc); acc = NULL; }
	    return;
	}
	
	/* allocate state if needed */
	if (R == NULL) {
	    /* mexPrintf("allocating R, Xp at %d\n", P+1); */
	    R = (float *)malloc((P+1)*sizeof(float));
	    for (i = 0; i <= P; ++i)  { R[i] = 0.0; }
	    R[0] = .001;
	    Xo = (float *)malloc((P+1)*sizeof(float));
	    for (i = 0; i < P; ++i)  { Xo[i] = 0.0; }
	    ai = (float *)malloc((P+1)*sizeof(float));
	    acc = (float *)malloc((P+1)*sizeof(float));
	}

	/* create output matrix */
	Ymatrix = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
	pY = mxGetData(Ymatrix);

	/* calculate autocorrelation of current block */
//	for (i = 0; i <= P; ++i) {
//	    double *x1 = pX;
//	    double *x2 = pX;
//	    float acc = 0;
//	    for (j = 0; j < N; ++j) {
//		if (j >= i) {
//		    acc += *x1 * *x2++;
//		}
//		++x1;
//	    }
//	    R[i] += alpha*(acc - R[i]);
//	}
	{
	    double *x1 = pX;
	    //float *acc = (float *)malloc((P+1)*sizeof(float));
	    for (i = 0; i<=P; ++i) { acc[i] = 0; }
	    for (j = 0; j < N; ++j) {
		double *x2 = x1;
		double x = *x1;
		float *pacc = acc;
		int lim = P;
		if (lim > j) lim = j;
		for (i = 0; i <= lim; ++i) {
		    *pacc++ += x * *x2--;
		}
		++x1;
	    }
	    /* smoothed update */
	    for (i = 0; i <=P; ++i) {
		R[i] += alpha*(acc[i] - R[i]);
	    }
	    //free(acc);
	}


	/* calculate new filter coefficients */
	/* Durbin's recursion, per p. 411 of Rabiner & Schafer 1978 */

	E = R[0];
	for (i = 1; i <= P; ++i) {
	    float sumalphaR = 0;
	    for (j = 1; j < i; ++j) {
		sumalphaR += ai[j]*R[i-j];
	    }
	    ki = (R[i] - sumalphaR)/E;
	    ai[i] = ki;
	    for (j = 1; j <= i/2; ++j) {
		float aj = ai[j];
		float aimj = ai[i-j];
		ai[j] = aj - ki*aimj;
		ai[i-j] = aimj - ki*aj;
	    }
	    E = (1-ki*ki)*E;
	}

	/* debug */
	if (0) {
	    for(i = 0; i <=P; ++i) {
		mexPrintf("%.5f ", R[i]);
	    }	
	    mexPrintf("\n");
	    for(i = 0; i <=P; ++i) {
		mexPrintf("%.5f ", ai[i]);
	    }
	}


	/* caculate new output */
	{
	    double *X = pX;
	    double *Xp = pX-1;
	    double *Y = pY;
	    float *pXo = Xo;
	    /* break loop over outputs into two blocks */
	    /* so that code that has to look back into history can be 
	       skipped for the bulk of the samples */
	    for (i = 0; i < P; ++i) {
		float acc = *X++;
		float *pai = ai+1;
		pXo = Xo + P;
		for (j = 1; j <= i; ++j) {
		    acc -= *pai++ * *Xp--;
		}
		for (j = i+1; j <= P; ++j) {
		    acc -= *pai++ * *(--pXo);
		}
		*Y++ = acc;
		Xp = X-1;
	    }
	    /* this is the bulk of the calculation - empirically fastest */
	    for (i = P; i < N; ++i) {
		//float acc = *X++;
		float acc = pX[i];
		//float *pai = ai+1;
		for (j = 1; j <= P; ++j) {
		    acc -= ai[j] * pX[i-j];  // fastest on mac !?
		    //acc -= *pai++ * *Xp--;
		    //acc -= *pai++ * pX[i-j];
		    //acc -= ai[j] * *Xp--;
		}
		*Y++ = acc;
		//pY[i] = acc;
		//Xp = X-1;
	    }
	}
	
	/* save last few frames of input */
	for (i = 0; i <= P; ++i) {
	    Xo[i] = pX[N-1-P+i];
	}

	plhs[0] = Ymatrix;
    }
}


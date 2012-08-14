/*
 * subband_analysis_helper.c
 * 
 * Matlab MEX implementation of Subband_Analysis_8.m
 *
 * Dan Ellis dpwe@ee.columbia.edu 2011-06-04
 */

#include "mex.h"
#include "string.h"
#include <math.h>

#ifndef MWSIZE_MIN
/* matlab R17 -ism? */
typedef unsigned int mwSize;
#endif

#define C_LEN 128
#define SUBBANDS 8

static float C[C_LEN] = {
 0.000000477,  0.000000954,  0.000001431,  0.000002384,  0.000003815,  0.000006199,  0.000009060,  0.000013828, 
 0.000019550,  0.000027657,  0.000037670,  0.000049591,  0.000062943,  0.000076771,  0.000090599,  0.000101566, 
-0.000108242, -0.000106812, -0.000095367, -0.000069618, -0.000027180,  0.000034332,  0.000116348,  0.000218868, 
 0.000339031,  0.000472546,  0.000611782,  0.000747204,  0.000866413,  0.000954151,  0.000994205,  0.000971317, 
-0.000868797, -0.000674248, -0.000378609,  0.000021458,  0.000522137,  0.001111031,  0.001766682,  0.002457142, 
 0.003141880,  0.003771782,  0.004290581,  0.004638195,  0.004752159,  0.004573822,  0.004049301,  0.003134727, 
-0.001800537, -0.000033379,  0.002161503,  0.004756451,  0.007703304,  0.010933399,  0.014358521,  0.017876148, 
 0.021372318,  0.024725437,  0.027815342,  0.030526638,  0.032754898,  0.034412861,  0.035435200,  0.035780907, 
-0.035435200, -0.034412861, -0.032754898, -0.030526638, -0.027815342, -0.024725437, -0.021372318, -0.017876148, 
-0.014358521, -0.010933399, -0.007703304, -0.004756451, -0.002161503,  0.000033379,  0.001800537,  0.003134727, 
-0.004049301, -0.004573822, -0.004752159, -0.004638195, -0.004290581, -0.003771782, -0.003141880, -0.002457142, 
-0.001766682, -0.001111031, -0.000522137, -0.000021458,  0.000378609,  0.000674248,  0.000868797,  0.000971317, 
-0.000994205, -0.000954151, -0.000866413, -0.000747204, -0.000611782, -0.000472546, -0.000339031, -0.000218868, 
-0.000116348, -0.000034332,  0.000027180,  0.000069618,  0.000095367,  0.000106812,  0.000108242,  0.000101566, 
-0.000090599, -0.000076771, -0.000062943, -0.000049591, -0.000037670, -0.000027657, -0.000019550, -0.000013828, 
-0.000009060, -0.000006199, -0.000003815, -0.000002384, -0.000001431, -0.000000954, -0.000000477, 0};

static float *Mr = NULL;
static float *Mi = NULL;

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int 	i,j,k;

    if (nrhs < 1){
	mexPrintf("U = subband_analysis_helper(X)\n");
	mexPrintf("           Perform 8-band subband analysis.\n");
	return;
    }

    if (nlhs > 0){
	mxArray  *Umatrix;
	double *pX, *pUr, *pUi;
	float *pY, *pZ;
	int t, i, j, k, xlen, hop, ucols;


	// Calculate the analysis filter bank coefficients
	if (Mr == NULL) {
#define M_ROWS 8
#define M_COLS 16
	    // mexPrintf("allocating Mr, Mi at %d points\n", M_ROWS*M_COLS);
	    Mr = (float *)malloc(M_ROWS*M_COLS*sizeof(float));
	    Mi = (float *)malloc(M_ROWS*M_COLS*sizeof(float));
	    for (i = 0; i < M_ROWS; ++i) {
		for (k = 0; k < M_COLS; ++k) {
		    Mr[i + M_ROWS*k] = cos((2*i + 1)*(k-4)*(M_PI/16.0));
		    Mi[i + M_ROWS*k] = sin((2*i + 1)*(k-4)*(M_PI/16.0));
		}
	    }
	}

	pZ = (float *)malloc(C_LEN*sizeof(float));
	pY = (float *)malloc(M_COLS*sizeof(float));

	xlen = mxGetM(prhs[0])*mxGetN(prhs[0]);
	pX = (double *)mxGetData(prhs[0]);

	if (xlen == 0) {
	    free(Mr);
	    free(Mi);
	    Mr = NULL;
	    Mi = NULL;
	}

	
 	/* mexPrintf("rows=%d cols=%d type=%s\n"  rows  cols  mxGetClassName(prhs[0])); */

	ucols = (xlen - C_LEN + 1)/SUBBANDS;

	Umatrix = mxCreateNumericMatrix(SUBBANDS, ucols, mxDOUBLE_CLASS, mxCOMPLEX);
	pUr = (double *)mxGetPr(Umatrix);
	pUi = (double *)mxGetPi(Umatrix);

	//mexPrintf("Alocated U with %d x %d\n", SUBBANDS, ucols);

	for (t = 0; t < ucols; ++t) {

	    //mexPrintf("t = %d\n", t);
	
	    for (i = 0; i < C_LEN; ++i) {
		pZ[i] = pX[ t*SUBBANDS + i] * C[i];
	    }

	    //mexPrintf("i1 = %d\n", i);
	    for (i = 0; i < M_COLS; ++i) {
		pY[i] = pZ[i];
	    }
	    //mexPrintf("i2 = %d\n", i);
	    for (i = 0; i < M_COLS; ++i) {
		for (j = 1; j < M_ROWS; ++j) {
		    pY[i] += pZ[i + M_COLS*j];
		}
	    }

	    //mexPrintf("i3 = %d\n", i);
	    for (i = 0; i < M_ROWS; ++i) {
		for (j = 0; j < M_COLS; ++j) {
		    pUr[(i + SUBBANDS*t)] += Mr[i + M_ROWS*j] * pY[j];
		    pUi[(i + SUBBANDS*t)] -= Mi[i + M_ROWS*j] * pY[j];
		}
	    }
	    //mexPrintf("i4 = %d\n", i);
	}

	free(pY);
	free(pZ);

	plhs[0] = Umatrix;
    }
}

#ifdef TEST
#endif /* TEST */

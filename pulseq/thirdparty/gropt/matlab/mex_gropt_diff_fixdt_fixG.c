#include <mex.h>
#include <stdio.h> 
#include <string.h>
#include "matrix.h"
#include "optimize_kernel.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *temp;
    double *G;
    int N;
    double *ddebug;
    const mwSize *dims;
    
    int verbose = 0;
    
    double gmax = mxGetScalar(prhs[0]);
    double smax = mxGetScalar(prhs[1]);
    
    int N_moments;
    double *moment_params;
    moment_params = mxGetPr(prhs[2]);
    dims = mxGetDimensions(prhs[2]);
    N_moments = dims[1];
    
    double TE = mxGetScalar(prhs[3]);
    double T_readout = mxGetScalar(prhs[4]);
    double T_90 = mxGetScalar(prhs[5]);
    double T_180 = mxGetScalar(prhs[6]);
    double dt = mxGetScalar(prhs[7]);
    int diffmode = mxGetScalar(prhs[8]);
    
    double dt_out;
    dt_out = mxGetScalar(prhs[9]);
    
    double pns_thresh;
    pns_thresh = mxGetScalar(prhs[10]);
    
    int N_eddy;
    double *eddy_params;
    eddy_params = mxGetPr(prhs[11]);
    dims = mxGetDimensions(prhs[11]);
    N_eddy = dims[1];
    
    double slew_reg = mxGetScalar(prhs[12]);
    
    int N_gfix;
    double *gfix;
    gfix = mxGetPr(prhs[13]);
    dims = mxGetDimensions(prhs[13]);
    N_gfix = dims[1];
    
    run_kernel_diff_fixeddt_fixG(&G, &N, &ddebug, verbose, dt, gmax, smax, TE, N_moments, moment_params, pns_thresh, 
            T_readout, T_90, T_180, diffmode, dt_out, N_eddy, eddy_params, -1.0, N_gfix, gfix, slew_reg);
    

    plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
    // Why doesnt this work?: mxSetPr(plhs[0], G);
    // Instead do this:
    temp = mxGetPr(plhs[0]);
    for (int i = 0; i < N; i++) {
        temp[i] = G[i];
    }
    
    plhs[1] = mxCreateDoubleScalar(ddebug[14]);
    
    // Need some frees here, m_params, eddy_params, G
        

}
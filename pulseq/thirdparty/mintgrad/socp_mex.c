/*
 *   socp_mex.c
 *
 *   Second-Order Cone Programming
 *   Interface to Matlab, compile with cmex
 *
 *   mlobo@isl.stanford.edu -- 96/97
 *   Modified by Nam Gyun Lee, namgyunl@usc.edu -- 01/02/2023
 */

#include <stdio.h>
#include <stddef.h>

#include "socp.h"
#include "mex.h"
#include "blas.h"

#include "socp.c"



/*
 *  [x, info, z, time, hist] =
 *       socp_mex(f, A, b, N, x, z, abs_tol, rel_tol, target, max_iter, Nu, out_mode)
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    ptrdiff_t n, m;
    ptrdiff_t L, *N;
    int iter, out_mode;
    double *f, *A, *b, *x, *z;
    double *x_opt, *z_opt;
    double abs_tol, rel_tol, target, Nu;
    double *hist;
    int mhist, nhist;
    int nptr;
    ptrdiff_t ndbl, nint;
    double *dblwork;
    ptrdiff_t *intwork;
    int i, info_time, info_exit;
    ptrdiff_t j, info_socp;
    double *info_ptr;
    static int firstcall = 1;
    double *time;

    ptrdiff_t int1 = 1;

    if (firstcall)
    {
        mexPrintf("\nThis is the beta version of SOCP.,\n");
        mexPrintf("COPYRIGHT (c) 1997, Miguel Lobo, Lieven Vandenberge, Stephen Boyd.\n\n");
        firstcall = 0;
    }

    /* check number of arguments */
    if (nrhs == 12)
    {
        if (nlhs != 5)
        {
            mexErrMsgTxt("Five output arguments required.\n");
        }
    } else
    {
        mexErrMsgTxt("Twelve input arguments required.\n");
    }

    /* get dimensions and ptr. to A */
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    A = mxGetPr(prhs[1]);

    /* check that dimensions of f, b, x and z match with A and get ptrs. */
    if (mxGetN(prhs[0]) != 1)
    {
        mexErrMsgTxt("1st input argument must be a column vector.");
    }
    if (mxGetM(prhs[0]) != n)
    {
        mexErrMsgTxt("f and A do not agree in size.");
    }
    f = mxGetPr(prhs[0]);
    
    if (mxGetN(prhs[2]) != 1)
    {
        mexErrMsgTxt("3rd input argument must be a column vector.");
    }
    if (mxGetM(prhs[2]) != m)
    {
        mexErrMsgTxt("A and b do not agree in size.");
    }
    b = mxGetPr(prhs[2]);

    if (mxGetN(prhs[4]) != 1)
    {
        mexErrMsgTxt("5th input argument must be a column vector.");
    }
    if (mxGetM(prhs[4]) != n)
    {
        mexErrMsgTxt("A and x do not agree in size.");
    }
    x = mxGetPr(prhs[4]);

    if (mxGetN(prhs[5]) != 1)
    {
        mexErrMsgTxt("6th input argument must be a column vector.");
    }
    if (mxGetM(prhs[5]) != m)
    {
        mexErrMsgTxt("A and z do not agree in size.");
    }
    z = mxGetPr(prhs[5]);

    /* get number of cones (=L) and size of each (=N[i]) */
    if (MIN(mxGetM(prhs[3]), mxGetN(prhs[3])) != 1)
    {
        mexErrMsgTxt("4th input argument must be a vector.\n");
    }
    L = MAX(mxGetM(prhs[3]), mxGetN(prhs[3]));
    N = mxCalloc(L, sizeof(ptrdiff_t)); /* new N of type ptrdiff_t */
    for (i = 0, j = 0; i < L; ++i)
    {
        N[i] = (ptrdiff_t) mxGetPr(prhs[3])[i];
        j += N[i];
    }
    if (j != m)
    {
        mexErrMsgTxt("Sum of elements of N must equal number of lines in A.\n");
    }

    /* get stopping criteria */
    abs_tol = mxGetScalar(prhs[6]);
    rel_tol = mxGetScalar(prhs[7]);
    target = mxGetScalar(prhs[8]);
    iter = (int) mxGetScalar(prhs[9]);

    /* get gap reduction vs. centering factor */
    Nu = mxGetScalar(prhs[10]);
    if (Nu < 0)
    {
        mexErrMsgTxt("Parameter Nu must be non-negative");
    }
    if (Nu <= 1)
    {
        fprintf(stdout, "SOCP warning: Nu > 1 recommended.\n");
    }

    /* output mode */
    out_mode = (int) mxGetScalar(prhs[11]);
    if (out_mode < 0 || out_mode > 2)
    {
        mexErrMsgTxt("Illegal value for out_mode; must be 0, 1 or 2.");
    }

    /* allocate workspace for socp() */
    socp_getwork(L, N, n, iter, out_mode, &mhist, &nhist, &ndbl, &nint);

    dblwork = (double *) mxCalloc(ndbl, sizeof(double));
    intwork = (ptrdiff_t *) mxCalloc(nint, sizeof(ptrdiff_t));

    /* prepare output arguments */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x_opt = mxGetPr(plhs[0]);
    dcopy_(&n, x, &int1, x_opt, &int1);  /* copy x to x_opt */

    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    info_ptr = mxGetPr(plhs[1]);

    plhs[2] = mxCreateDoubleMatrix(m, 1, mxREAL);
    z_opt = mxGetPr(plhs[2]);
    dcopy_(&m, z, &int1, z_opt, &int1);  /* copy z to z_opt */

    plhs[3] = mxCreateDoubleMatrix(mhist, nhist, mxREAL);
    hist = mxGetPr(plhs[3]);

    plhs[4] = mxCreateDoubleMatrix(1, 3, mxREAL);  /* time stats: utime, stime, iters */
    time = mxGetPr(plhs[4]);

    /* call socp() */
    info_socp = socp(L, N, n, f, A, b, x_opt, z_opt,
                     abs_tol, rel_tol, target, &iter, Nu,
                     &info_exit, out_mode, hist, dblwork, intwork);

    time[2] = iter;

    *info_ptr = (double) info_exit;

    /* rescale hist */

    /* note: hist will padded with zeros after call to socp; the following */
    /* line truncates it */
    /* (if hist is used to output other info, this may need to be removed) */
    mxSetN(plhs[3], iter + 1);

    /* free allocated memory */
    mxFree(N);
    mxFree(dblwork);
    mxFree(intwork);

    /* error handling */
    if (info_socp)
    {
        mexPrintf("info from dgelss = %d \n", info_socp);
        mexErrMsgTxt("Error in SOCP, call to LAPACK failed.\n");
    }
}

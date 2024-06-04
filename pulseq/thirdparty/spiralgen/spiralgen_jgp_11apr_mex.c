/*=======================================================================*/
/*                                                                       */
/*  SOURCE_FILE:    SPIRALGEN_MEX.C                                      */
/*                                                                       */
/*  Copyright 2016: Namgyun Lee                                          */
/*  Author: Namgyun Lee                                                  */
/*  USC (University of Southern California)                              */
/*  KBSI (Korea Basic Science Institute)                                 */
/*  namgyunl@kbsi.re.kr, ggang56@gmail.com (preferred)                   */
/*=======================================================================*/
/*=======================================================================*/
/*  I N C L U D E S                                                      */
/*=======================================================================*/
/*-----------------------------------------------------------------------*/
/*  Common                                                               */
/*-----------------------------------------------------------------------*/
#define _USE_MATH_DEFINES /* To use MATH Constants */
#include <math.h>

/*-----------------------------------------------------------------------*/
/*  MATLAB                                                               */
/*-----------------------------------------------------------------------*/
#include "mex.h"
#include "spiralgen_jgp_11apr.c"

/*=======================================================================*/
/*  G L O B A L   R E F E R E N C E S                                    */
/*=======================================================================*/
/*=======================================================================*/
/*  G L O B A L   D E F I N I T I O N S                                  */
/*=======================================================================*/

/*=======================================================================*/
/*  L O C A L   S Y M B O L   D E F I N I T I O N S                      */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   D A T A   D E F I N I T I O N S                          */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   F U N C T I O N   P R O T O T Y P E S                    */
/*=======================================================================*/
/*=======================================================================*/
/*  F U N C T I O N   P R O T O T Y P E S                                */
/*=======================================================================*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *spparams = NULL;

    double gamma;    /* typically 4257 [Hz/G] */
    double fov;      /* field of view [cm] */
    double res;      /* resolution [cm[ */
    double slewmax;  /* max slew rate [G/cm/sec] */
    double greadmax; /* max gradient amplitude (system or bw-limited) [G/cm] */
    double gsysmax;  /* total max gradient amplitude [G/cm] */
    double rast;     /* "raster time" between samples in sec, i.e. 0.000004 gives 4 usec gradient update */
    int    arms;     /* number of spiral interleaves */

    double us_0;
    double us_1;
    double us_r;
    int    ustype; /* rate of change in undersampling (see code below) */
                   /* SET us_0 = 1 or us_r = 1 to avoid variable sampling */

    /* some more variables determining the type of waveform */
    int gtype; /* 0 = calculate through readout, 1 = include grad rampdown, 2 = include rewinder to end at k=0 */
    int sptype; /* 0 = Archimedean, 1 = Fermat */

    size_t maxarray = 50000; /* maximum size of array -change if necessary */
    double *karray;
    double *garray;
    int karrlen;
    int garrlen;

    /*-------------------------------------------------------------------*/
    /* Check for proper number of arguments                              */
    /*-------------------------------------------------------------------*/
    if (nrhs != 1)
    {
        mexErrMsgTxt("1 input is required.");
    }
    if (nlhs != 4)
    {
        mexErrMsgTxt("4 outputs are required.");
    }

    /*-------------------------------------------------------------------*/
    /* Check the data type of inputs                                     */
    /*-------------------------------------------------------------------*/
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    {
        mexErrMsgTxt("The 1st input, spparams, must be double and real.");
    }

    /*-------------------------------------------------------------------*/
    /* Parse the input parameters                                        */
    /*-------------------------------------------------------------------*/
    spparams = (double *) mxGetData(prhs[0]);
    gamma    = spparams[0];
    fov      = spparams[1];
    res      = spparams[2];
    slewmax  = spparams[3];
    greadmax = spparams[4];
    gsysmax  = spparams[5];
    rast     = spparams[6];
    arms     = (int) spparams[7];
    us_0     = spparams[8];
    us_1     = spparams[9];
    us_r     = spparams[10];
    ustype   = (int) spparams[11];
    gtype    = (int) spparams[12];
    sptype   = (int) spparams[13];

    printf("gamma    = %f\n", gamma);
    printf("fov      = %f\n", fov);
    printf("res      = %f\n", res);
    printf("slewmax  = %f\n", slewmax);
    printf("greadmax = %f\n", greadmax);
    printf("gsysmax  = %f\n", gsysmax);
    printf("rast     = %f\n", rast);
    printf("arms     = %d\n", arms);
    printf("us_0     = %f\n", us_0);
    printf("us_1     = %f\n", us_1);
    printf("us_r     = %f\n", us_r);
    printf("ustype   = %d\n", ustype);
    printf("gtype    = %d\n", gtype);
    printf("sptype   = %d\n", sptype);

    /*-------------------------------------------------------------------*/
    /* Create output MATLAB arrays for the return argument               */
    /*-------------------------------------------------------------------*/
    plhs[0] = mxCreateDoubleMatrix(2 * maxarray, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2 * maxarray, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    karray = (double *) mxGetData(plhs[0]);
    garray = (double *) mxGetData(plhs[1]);

    /*-------------------------------------------------------------------*/
    /* Call spiralgen_jgp_11apr.c                                        */
    /*-------------------------------------------------------------------*/
    spiralgen_jgp_11apr(spparams, karray, garray, &karrlen, &garrlen);

    printf("karrlen = %d\n", karrlen);
    printf("garrlen = %d\n", garrlen);

    /*-------------------------------------------------------------------*/
    /* Assign outputs                                                    */
    /*-------------------------------------------------------------------*/
    *(mxGetPr(plhs[2])) = karrlen;
    *(mxGetPr(plhs[3])) = garrlen;
}

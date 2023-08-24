#include <System/SizedIntegers.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>


#include "AlgAnisotropicGaussHybridCubic.h"




inline double interpolateCubicly(double a, double b, double c, double d, double y0, double y1, double y20, double y21)
{
    return a*y0 + b*y1 + c*y20 + d*y21;
}



AlgAnisotropicGaussHybridCubic::AlgAnisotropicGaussHybridCubic(double *input, int sizeX, int sizeY)
    : AlgAnisotropicGaussBase(input, sizeX, sizeY)
{
}


AlgAnisotropicGaussHybridCubic::~AlgAnisotropicGaussHybridCubic()
{
}



void AlgAnisotropicGaussHybridCubic::runFilter(double *output, double sigma_u, double sigma_v, double theta)
{
    // correct the angle
    theta = 180.0-theta;
    computeDeviationsAndIntercept(sigma_u, sigma_v, theta, false);

    // filter wrt x-axis
    filterXAxis(m_input, output, output, m_sigma_x);
    // filter wrt t-line (or y-axis if no interpolation needed)
    if(theta == 180 || theta == 90)
    {
        filterYAxis(output, output, output, m_sigma_phi);
    }
    else
    {
        filterTLine(output, output, output);
    }
}



void AlgAnisotropicGaussHybridCubic::calculateY2(int x, double *y, double *y2, double *u)
{
    double p;
    int sizeX = m_sizeX;
    y2[0] = 0.0;
    u[0] = 0.0;
    for (int i = 1; i < sizeX - 1; i++)
    {
        p = y2[i - 1]/2.0 + 2.0;
        y2[i] = -1.0/(2.0*p);
        u[i] = y[i + 1 + x*sizeX] + y[i - 1 + x*sizeX] - 2.0*y[i + x*sizeX];
        u[i] = (3.0*u[i] - u[i - 1]/2.0)/p;
    }
    y2[sizeX - 1] = 0.0;
    u[sizeX - 1] = 0.0;
    for (int k = sizeX-2; k >= 0; k--)
    {
        y2[k] = y2[k]*y2[k + 1] + u[k];
    }
}



// --------------- filter t-line ------------------------------------------

void AlgAnisotropicGaussHybridCubic::filterTLine(double *input, double *w, double *output)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;

    double tan_phi = m_tan_phi;
    double yMu;
    // setup
    filterLines wLines;
    wLines.line_n = new double [sizeX + int(sizeY/fabs(tan_phi)) + 2];
    wLines.line_n1 = new double [sizeX + int(sizeY/fabs(tan_phi)) + 2];
    wLines.line_n2 = new double [sizeX + int(sizeY/fabs(tan_phi)) + 2];
    wLines.line_n3 = new double [sizeX + int(sizeY/fabs(tan_phi)) + 2];

    double *lastLineOfInput = new double[sizeX];
    for (int x = 0; x < sizeX; x++)
    {
        lastLineOfInput[x] = input[(sizeY - 1)*sizeX + x];
    }

    //needed for interpolation
    double* y2 = new double [sizeX];
    double* u2 = new double [sizeX];
    double a, b, c, d;
    double prevEntry, currentEntry;

    // if tan(phi) is negative, the t-lines run out of the image on the left side
    // --> make sure that memory is allocated there
    if(tan_phi < 0)
    {
        wLines.line_n = &(wLines.line_n)[int(sizeY/fabs(tan_phi) + 1)];
        wLines.line_n1 = &(wLines.line_n1)[int(sizeY/fabs(tan_phi) + 1)];
        wLines.line_n2 = &(wLines.line_n2)[int(sizeY/fabs(tan_phi) + 1)];
        wLines.line_n3 = &(wLines.line_n3)[int(sizeY/fabs(tan_phi) + 1)];
    }

    coefficients coeffs;
    double q = calculateQ(m_sigma_phi);
    calculateBCoefficients(q, &coeffs);
    calculateM(&coeffs);

    //interpolation step before filtering
    for (int x = 1; x < sizeY; x++)
    {
        yMu = double(x)/tan_phi;

        calculateY2(x, input, y2, u2);
        //interpolation parameters
        a = yMu - floor(yMu);
        b = 1-a;
        c = (a*a*a - a)/6.0;
        d = (b*b*b - b)/6.0;

        currentEntry = input[x*sizeX];

        for (int y = 1; y < sizeX;  y++)
        {
            prevEntry = currentEntry;
            currentEntry = input[y + x*sizeX];
            w[y + x*sizeX] = interpolateCubicly(a,b,c,d,prevEntry, currentEntry, y2[y - 1], y2[y]);
        }

    }
    filterYoungForwardT(w, w, wLines, coeffs);

    // these swaps within filter_young_forward_t don't happen outside...
    // in order to find the current w_ni, swap here too
    for(int i = 0; i < (sizeY - 1)%4; i++)
    {
        double *tmpSwap = wLines.line_n3; // remember where we allocated that piece of memory, so we can reuse it in the next iteration
        wLines.line_n3 = wLines.line_n2;
        wLines.line_n2 = wLines.line_n1;
        wLines.line_n1 = wLines.line_n;
        wLines.line_n = tmpSwap;
    }

    filterYoungBackwardT(lastLineOfInput, w, w, wLines, coeffs);

    //interpolate back
    for (int x = sizeY - 2; x >= 0; x--)
    {
        yMu = double(x)/tan_phi;

        calculateY2(x, w, y2, u2);
        //interpolation parameters
        a = yMu - floor(yMu);
        b = 1-a;
        c = (a*a*a - a)/6.0;
        d = (b*b*b - b)/6.0;

        currentEntry = w[sizeX - 1 + x*sizeX];

        for (int y = sizeX - 2; y >= 0; y--)
        {
            prevEntry = currentEntry;
            currentEntry = w[y + x*sizeX];
            output[y + x*sizeX] = interpolateCubicly(a,b,c,d,prevEntry, currentEntry, y2[y + 1], y2[y]);
        }

    }


    // delete all allocated memory
    if(tan_phi < 0)
    {
        wLines.line_n = &(wLines.line_n)[ - int(sizeY/fabs(tan_phi) + 1)];
        wLines.line_n1 = &(wLines.line_n1)[ - int(sizeY/fabs(tan_phi) + 1)];
        wLines.line_n2 = &(wLines.line_n2)[ - int(sizeY/fabs(tan_phi) + 1)];
        wLines.line_n3 = &(wLines.line_n3)[ - int(sizeY/fabs(tan_phi) + 1)];
    }

    delete [] lastLineOfInput;
    delete [] wLines.line_n;
    delete [] wLines.line_n1;
    delete [] wLines.line_n2;
    delete [] wLines.line_n3;

    delete [] y2;
    delete [] u2;

}



// ---- filter forward direction t-line
// We iterate now over lines in x-direction.
// --> more memory, but we need to buffer anyway
// --> more efficient since we address neighboring values (?)
void AlgAnisotropicGaussHybridCubic::filterYoungForwardT(double *input, double* w, filterLines wLines, coefficients coeffs)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double tan_phi = m_tan_phi;
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};
    double sqrtb4 = sqrt(b[4]);
    double yMu;
    int yMuDiscrete, yMuDiscretePrev = 0;
    double* w_n = wLines.line_n;
    double* w_n1 = wLines.line_n1;
    double* w_n2 = wLines.line_n2;
    double* w_n3 = wLines.line_n3;

    // prepare values for first iteration
    for (int y = 0; y < sizeX; y++)
    {
        //w[y] = input[y];// initial values //we already did that in the step before
        w_n1[y] = w[y];
        w_n2[y] = w[y];
        w_n3[y] = w[y];

    }

    for (int x = 1; x < sizeY; x++)
    {
        yMu = double(x)/tan_phi;
        yMuDiscrete = int(yMu);
        // handle that t-lines will move, so we need to fill up empty indices
        if(tan_phi > 0)
        {
            for (int y = 0; y < yMuDiscrete - yMuDiscretePrev; y++)
            {
                w_n1[sizeX + yMuDiscretePrev + y] = w_n1[sizeX + yMuDiscretePrev - 1];
                w_n2[sizeX + yMuDiscretePrev + y] = w_n2[sizeX + yMuDiscretePrev - 1];
                w_n3[sizeX + yMuDiscretePrev + y] = w_n3[sizeX + yMuDiscretePrev - 1];
            }
        }
        else
        {
            for (int y = 0; y < abs(yMuDiscrete) - abs(yMuDiscretePrev); y++)
            {
                w_n1[yMuDiscretePrev - y - 1] = w_n1[yMuDiscretePrev];
                w_n2[yMuDiscretePrev - y - 1] = w_n2[yMuDiscretePrev];
                w_n3[yMuDiscretePrev - y - 1] = w_n3[yMuDiscretePrev];
            }
        }

        w_n[yMuDiscrete] = b[1]*w_n1[yMuDiscrete] + b[2]*w_n2[yMuDiscrete] + b[3]*w_n3[yMuDiscrete];
        w_n[yMuDiscrete] = sqrtb4*input[x*sizeX] - w_n[yMuDiscrete];
        w[x*sizeX] = w_n[yMuDiscrete];

        for (int y = 1; y < sizeX; y++)
        {
            // step
            w_n[y + yMuDiscrete] = b[1]*w_n1[y + yMuDiscrete] + b[2]*w_n2[y + yMuDiscrete] + b[3]*w_n3[y + yMuDiscrete];
            w_n[y + yMuDiscrete] = sqrtb4*input[y + x*sizeX] - w_n[y + yMuDiscrete];
            w[y + x*sizeX] = w_n[y + yMuDiscrete];
        }

        // swap history
        double *tmpSwap = w_n3; // remember where we allocated that piece of memory, so we can reuse it in the next iteration
        w_n3 = w_n2;
        w_n2 = w_n1;
        w_n1 = w_n;
        w_n = tmpSwap;

        yMuDiscretePrev = yMuDiscrete;

    }
}


// ---- filter backward direction t-line
void AlgAnisotropicGaussHybridCubic::filterYoungBackwardT(double* lastLineOfInput, double* w, double* output, filterLines outputLines, coefficients coeffs)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double tan_phi = m_tan_phi;
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};
    double sqrtb4 = sqrt(b[4]);
    double yMu = double(sizeY - 1)/m_tan_phi;
    int yMuDiscrete, yMuDiscretePrev = int(yMu);
    double *tmpSwap;
    double* output_n = outputLines.line_n;
    double* output_n1 = outputLines.line_n1;
    double* output_n2 = outputLines.line_n2;
    double* output_n3 = outputLines.line_n3;

    // prepare values for first iteration (incl Triggs)
    filterYoungBackwardInitialValuesT(lastLineOfInput, output, outputLines, coeffs);

    // iteration (incl t-line handling)
    for (int x = sizeY - 2; x >= 0; x--)
    {
        yMu = double(x)/tan_phi;
        yMuDiscrete = int(yMu);
        // handle that t-lines will move, so we need to fill up empty indices
        if(tan_phi > 0)
        {
            for (int y = 0; y < yMuDiscretePrev - yMuDiscrete; y++)
            {
                output_n1[yMuDiscretePrev - y - 1] = output_n1[yMuDiscretePrev];
                output_n2[yMuDiscretePrev - y - 1] = output_n2[yMuDiscretePrev];
                output_n3[yMuDiscretePrev - y - 1] = output_n3[yMuDiscretePrev];
            }
        }
        else
        {
            for(int y = 0; y < abs(yMuDiscretePrev) - abs(yMuDiscrete); y++)
            {
                output_n1[sizeX + yMuDiscretePrev + y] = output_n1[sizeX  - 1 + yMuDiscretePrev];
                output_n2[sizeX + yMuDiscretePrev + y] = output_n2[sizeX  - 1 + yMuDiscretePrev];
                output_n3[sizeX + yMuDiscretePrev + y] = output_n3[sizeX  - 1 + yMuDiscretePrev];
            }
        }
        output_n[sizeX - 1 + yMuDiscrete] = b[1]*output_n1[sizeX - 1 + yMuDiscrete] + b[2]*output_n2[sizeX - 1 + yMuDiscrete] + b[3]*output_n3[sizeX - 1 + yMuDiscrete];
        output_n[sizeX - 1 + yMuDiscrete] = sqrtb4*w[sizeX - 1 + x*sizeX] - output_n[sizeX - 1 + yMuDiscrete];
        output[sizeX - 1 + x*sizeX] = output_n[sizeX - 1 + yMuDiscrete];

        for (int y = sizeX - 2; y >= 0; y--)
        {
            // step
            output_n[y + yMuDiscrete] = b[1]*output_n1[y + yMuDiscrete] + b[2]*output_n2[y + yMuDiscrete] + b[3]*output_n3[y + yMuDiscrete];
            output_n[y + yMuDiscrete] = sqrtb4*w[y + x*sizeX] - output_n[y + yMuDiscrete];
            output[y + x*sizeX] = output_n[y + yMuDiscrete];
        }

        // swap history
        tmpSwap = output_n3; // remember where we allocated that piece of memory, so we can reuse it in the next iteration
        output_n3 = output_n2;
        output_n2 = output_n1;
        output_n1 = output_n;
        output_n = tmpSwap;

        yMuDiscretePrev = yMuDiscrete;
    }

}



void AlgAnisotropicGaussHybridCubic::filterYoungBackwardInitialValuesT(double *lastLineOfInput, double* output, filterLines outputLines, coefficients coeffs)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double m11 = coeffs.M[0], m12 = coeffs.M[1], m13 = coeffs.M[2];
    double m21 = coeffs.M[3], m22 = coeffs.M[4], m23 = coeffs.M[5];
    double m31 = coeffs.M[6], m32 = coeffs.M[7], m33 = coeffs.M[8];
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};
    double sqrtb4 = sqrt(b[4]);

    double* output_n1 = outputLines.line_n1;
    double* output_n2 = outputLines.line_n2;
    double* output_n3 = outputLines.line_n3;

    double u_plus, v_plus;
    double u_t, u_t1, u_t2;


    // Triggs05

    // for better readability, name the variables such that they correspond
    // to the variable names in Triggs05
    double b1 = -b[1];
    double b2 = -b[2];
    double b3 = -b[3];


    double yMu = double(sizeY - 1.0)/m_tan_phi;

    for (int y = sizeX - 1; y >= 0; y--)
    {
        u_plus = lastLineOfInput[y]*sqrtb4/(1.0 - b1 - b2 - b3);
        v_plus = u_plus/(1.0 - b1 - b2 - b3);

        u_t = output_n1[y + int(yMu)] - u_plus;
        u_t1 = output_n2[y + int(yMu)] - u_plus;
        u_t2 = output_n3[y + int(yMu)] - u_plus;

        output_n1[y + int(yMu)] = (m11*u_t + m12*u_t1 + m13*u_t2 + v_plus)*sqrtb4;
        output[(sizeY - 1)*sizeX + y] = output_n1[y + int(yMu)];
        output_n2[y + int(yMu)] = (m21*u_t + m22*u_t1 + m23*u_t2 + v_plus)*sqrtb4;
        output_n3[y + int(yMu)] = (m31*u_t + m32*u_t1 + m33*u_t2 + v_plus)*sqrtb4;
    }

}


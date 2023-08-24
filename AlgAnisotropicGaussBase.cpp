#include <System/SizedIntegers.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>


#include "AlgAnisotropicGaussBase.h"




AlgAnisotropicGaussBase::AlgAnisotropicGaussBase(double *input, int sizeX, int sizeY)
    : m_sizeX(sizeX)
    , m_sizeY(sizeY)
    , m_sigma_x(0.0)
    , m_sigma_phi(0.0)
    , m_tan_phi(0.0)
    , m_input(input)
{
    m_initialValuesBackwards = new double[3];
}


AlgAnisotropicGaussBase::~AlgAnisotropicGaussBase()
{
    delete [] m_initialValuesBackwards;
}




// ---- calculate coefficients

void AlgAnisotropicGaussBase::computeDeviationsAndIntercept(double sigma_u, double sigma_v, double theta, bool flipXY)
{
    if (!flipXY)
    {
        double theta_rad = theta*M_PI/180;

        double sigma_u2 = sigma_u*sigma_u;
        double sigma_v2 = sigma_v*sigma_v;
        double squares_and_trigs = sigma_v2*cos(theta_rad)*cos(theta_rad);
        squares_and_trigs += sigma_u2*sin(theta_rad)*sin(theta_rad);

        // standard deviation x-axis: formula (9)
        m_sigma_x = sigma_u*sigma_v/sqrt(squares_and_trigs);
        // standard deviation t-line: formula (10)
        m_sigma_phi = sqrt(squares_and_trigs); // omit /sin(phi) because it vanishes in (4)
        // intercept for t-line: formula (11)
        m_tan_phi = squares_and_trigs/((sigma_u2 - sigma_v2)*cos(theta_rad)*sin(theta_rad));
    }
    else
    {
        double theta_rad = theta*M_PI/180;

        double sigma_u2 = sigma_u*sigma_u;
        double sigma_v2 = sigma_v*sigma_v;
        double squares_and_trigs = sigma_v2*cos(M_PI/2.0 + theta_rad)*cos(M_PI/2.0 + theta_rad);
        squares_and_trigs += sigma_u2*sin(M_PI/2.0 + theta_rad)*sin(M_PI/2.0 + theta_rad);

        // standard deviation x-axis: formula (9)
        m_sigma_x = sigma_u*sigma_v/sqrt(squares_and_trigs);
        // standard deviation t-line: formula (10)
        m_sigma_phi = sqrt(squares_and_trigs); // omit /sin(phi) because it vanishes in (4)
        // intercept for t-line: formula (11)
        m_tan_phi = squares_and_trigs/((sigma_u2 - sigma_v2)*cos(M_PI/2.0 +theta_rad)*sin(M_PI/2.0 +theta_rad));
        m_tan_phi = tan(M_PI - atan(m_tan_phi));
    }

}


double AlgAnisotropicGaussBase::calculateQ(double sigma)
{

    // These values were determined in Young02 for the Gabor filter.
    // Since the Gabor filter and Gauss filter relate through modulation, the values
    // should hold for the Gauss filter too. Anyway, I do not see any reason why
    // these values should be better than those computed in Young95
    if (sigma >= 3.556)
    {
        return 2.5091 + 0.9804*(sigma - 3.556);
    }
    else if (sigma >= 0.5)
    {
        return -0.2568 + 0.5784*sigma + 0.0561*sigma*sigma;
    }
    else
    {
        std::cout << "WARNING: The Gaussian envelope is too narrow!" << std::endl;
        return sigma;
    }

}


void AlgAnisotropicGaussBase::calculateBCoefficients(double q, coefficients* coeffs)
{
    // Young02 equation (6)
    double m_0 = 1.16680;
    double m_1 = 1.10783;
    double m_2 = 1.40586;

    // Young02 equation (10)
    coeffs->b[0] = 1.0;

    double scale = m_1*m_1 + m_2*m_2 + 2*m_1*q + q*q;
    scale = (m_0+q)*scale;

    coeffs->b[1] = 2*m_0*m_1 + m_1*m_1 + m_2*m_2 + (2*m_0 + 4*m_1)*q + 3*q*q;
    coeffs->b[1] = -q*coeffs->b[1]/scale;

    coeffs->b[2] = m_0 + 2*m_1 + 3*q;
    coeffs->b[2] = q*q*coeffs->b[2]/scale;

    coeffs->b[3] = -q*q*q/scale;
    coeffs->b[4] = m_0*(m_1*m_1 + m_2*m_2)/scale;
    coeffs->b[4] = coeffs->b[4]*coeffs->b[4];

}


/*
 * This computes the matrix M as defined in Triggs05. The matrix is necessary
 * to compute the initial values for the Gaussian filters in order to avoid distortions
 * due to discretisation.
 */
void AlgAnisotropicGaussBase::calculateM(coefficients* coeffs)
{

    // for better readability, name the variables such that they correspond
    // to the variable names in Triggs05
    double a1 = -coeffs->b[1];
    double a2 = -coeffs->b[2];
    double a3 = -coeffs->b[3];

    // matrix M (eq. 3) in Triggs05
    double scale = (1.0 + a1 - a2 + a3);
    scale *= (1.0 - a1 - a2 - a3);
    scale *= (1.0 + a2 + (a1 - a3)*a3);

    coeffs->M[0] = (-a3*a1 + 1.0 - a3*a3 - a2)/scale;
    coeffs->M[1] = (a3 + a1)*(a2 + a3*a1)/scale;
    coeffs->M[2] = a3*(a1 + a3*a2)/scale;

    coeffs->M[3] = (a1 + a3*a2)/scale;
    coeffs->M[4] = -(a2 - 1.0)*(a2 + a3*a1)/scale;
    coeffs->M[5] = -(a3*a1 + a3*a3 + a2 -1.0)*a3/scale;

    coeffs->M[6] = (a3*a1 + a2 + a1*a1 - a2*a2)/scale;
    coeffs->M[7] = (a1*a2 + a3*a2*a2 - a1*a3*a3 - a3*a3*a3 - a3*a2 + a3)/scale;
    coeffs->M[8] = a3*(a1 + a3*a2)/scale;
}



// ------- filter direction x-axis -------------------
void AlgAnisotropicGaussBase::filterXAxis(double *input, double *w, double *output, double sigma)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    // save input for later as it will be overwritten by filter_young_forward_x
    double *lastColumnOfInput = new double[sizeY];
    for (int y = 0; y < sizeY; y++)
    {
        lastColumnOfInput[y] = input[(sizeX - 1) + y*sizeX];
    }

    //double q = calculateQ(m_sigma_x); // an adapted standard deviation
    double q = calculateQ(sigma); // an adapted standard deviation
    coefficients coeffs;
    calculateBCoefficients(q, &coeffs); // coefficients for recursive equation
    calculateM(&coeffs); // necessary for computations of initial values

    filterYoungForwardX(input, w, coeffs);
    filterYoungBackwardX(input, lastColumnOfInput, w, output, coeffs);

    delete[] lastColumnOfInput;
}


// ---- filter forward direction x-axis
void AlgAnisotropicGaussBase::filterYoungForwardX(double *input, double* w, coefficients coeffs)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};
    double w_n, w_n1, w_n2, w_n3;

    for (int y = 0; y < sizeY; y++)
    {
        // initial values ("outside" of the image)
        w_n1 = input[0 + y*sizeX]/(1.0 + b[1] + b[2] + b[3]);
        w_n2 = w_n1;
        w_n3 = w_n1;
        //iterate over image pixels in x-dir
        for (int x = 0; x < sizeX; x++)
        {
            // step
            w_n = b[1]*w_n1 + b[2]*w_n2 + b[3]*w_n3;
            w_n = (input[x + y*sizeX] - w_n);

            w[x + y*sizeX] = w_n;
            w_n3 = w_n2;
            w_n2 = w_n1;
            w_n1 = w_n;
        }
    }
}


// ---- filter backward direction x-axis

void AlgAnisotropicGaussBase::filterYoungBackwardX(double* input,
                                                   double *lastColumnOfInput,
                                                   double* w, double* output,
                                                   coefficients coeffs)
{
    double result;
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};

    for (int y = 0; y < sizeY; y++)
    {
        // get initial values ("outside" of the image)
        filterYoungBackwardInitialValuesX(input, lastColumnOfInput, w, y, coeffs);
        double output_n3;
        double output_n2 = m_initialValuesBackwards[2];
        double output_n1 = m_initialValuesBackwards[1];
        double output_n = m_initialValuesBackwards[0];

        // prepare values for loop
        output[(sizeX - 1) + y*sizeX] = output_n;
        output_n3 = output_n2;
        output_n2 = output_n1;
        output_n1 = output_n;

        // iterate over image pixels in y-dir
        for (int x = sizeX - 2; x >= 0; x--)
        {
            // step
            result = b[1]*output_n1 + b[2]*output_n2 + b[3]*output_n3;
            output[x + y*sizeX] = (b[4]*w[x + y*sizeX] - result);

            // swap history
            output_n3 = output_n2;
            output_n2 = output_n1;
            output_n1 = output[x + y*sizeX];
        }
    }
}


void AlgAnisotropicGaussBase::filterYoungBackwardInitialValuesX(double *input,
                                                                double *lastColumnOfInput,
                                                                double* w,
                                                                int offset,
                                                                coefficients coeffs)
{
    int sizeX = m_sizeX;

    double m11 = coeffs.M[0], m12 = coeffs.M[1], m13 = coeffs.M[2];
    double m21 = coeffs.M[3], m22 = coeffs.M[4], m23 = coeffs.M[5];
    double m31 = coeffs.M[6], m32 = coeffs.M[7], m33 = coeffs.M[8];
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};


    // Triggs05

    // for better readability, name the variables such that they correspond
    // to the variable names in Triggs05
    double a1 = -b[1];
    double a2 = -b[2];
    double a3 = -b[3];
    double b1 = -b[1];
    double b2 = -b[2];
    double b3 = -b[3];

    // eq. 15
    double i_plus = lastColumnOfInput[offset];
    double u_plus = i_plus/(1.0 - a1 - a2 - a3);
    double v_plus = u_plus/(1.0 - b1 - b2 - b3);

    double u_t = w[(sizeX - 1) + offset*sizeX] - u_plus;
    double u_t1 = w[(sizeX - 2) + offset*sizeX] - u_plus;
    double u_t2 = w[(sizeX - 3) + offset*sizeX] - u_plus;

    m_initialValuesBackwards[0] = m11*u_t + m12*u_t1 + m13*u_t2 + v_plus;
    m_initialValuesBackwards[0] *= b[4];
    m_initialValuesBackwards[1] = m21*u_t + m22*u_t1 + m23*u_t2 + v_plus;
    m_initialValuesBackwards[1] *= b[4];
    m_initialValuesBackwards[2] = m31*u_t + m32*u_t1 + m33*u_t2 + v_plus;
    m_initialValuesBackwards[2] *= b[4];

}


// --------------- filter y-axis ----------------------------------------
void AlgAnisotropicGaussBase::filterYAxis(double *input, double *w, double *output, double sigma)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    filterLines wLines;
    wLines.line_n = new double [sizeX];
    wLines.line_n1 = new double [sizeX];
    wLines.line_n2 = new double [sizeX];
    wLines.line_n3 = new double [sizeX];

    // save input for later as it will be overwritten by filter_young_forward_y
    double *lastLineOfInput = new double[sizeX];
    for (int x = 0; x < sizeX; x++)
    {
        lastLineOfInput[x] = input[(sizeY - 1)*sizeX + x];
    }

    //double q = calculateQ(m_sigma_phi);
    double q = calculateQ(sigma);
    coefficients coeffs;
    calculateBCoefficients(q, &coeffs);
    calculateM(&coeffs);
    filterYoungForwardY(input, w, wLines, coeffs);

    // these swaps within filter_young_forward_y don't happen outside...
    // in order to find the current w_ni, swap here too
    for(int i = 0; i < (sizeY - 1)%4; i++)
    {
        double *tmpSwapp = wLines.line_n3; // remember where we allocated that piece of memory, so we can reuse it in the next iteration
        wLines.line_n3 = wLines.line_n2;
        wLines.line_n2 = wLines.line_n1;
        wLines.line_n1 = wLines.line_n;
        wLines.line_n = tmpSwapp;
    }

    filterYoungBackwardY(lastLineOfInput, w, output, wLines, coeffs);

    delete [] lastLineOfInput;
    delete [] wLines.line_n;
    delete [] wLines.line_n1;
    delete [] wLines.line_n2;
    delete [] wLines.line_n3;

}


// ---- filter forward direction y-axis
void AlgAnisotropicGaussBase::filterYoungForwardY(double *input, double* w, filterLines wLines, coefficients coeffs)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};
    double sqrtb4 = sqrt(b[4]);

    double* w_n = wLines.line_n;
    double* w_n1 = wLines.line_n1;
    double* w_n2 = wLines.line_n2;
    double* w_n3 = wLines.line_n3;

    // prepare values for first iteration
    for (int x = 0; x < sizeX; x++)
    {
        w[x] = input[x];// initial values
        w_n1[x] = w[x];
        w_n2[x] = w[x];
        w_n3[x] = w[x];

    }


    for (int y = 1; y < sizeY; y++)
    {

        w_n[0] = b[1]*w_n1[0] + b[2]*w_n2[0] +b[3]*w_n3[0];
        w_n[0] = sqrtb4*input[y*sizeX] - w_n[0];
        w[y*sizeX] = w_n[0];

        for (int x = 1; x < sizeX; x++)
        {
            // step
            w_n[x] = b[1]*w_n1[x] + b[2]*w_n2[x] + b[3]*w_n3[x];
            w_n[x] = sqrtb4*input[x + y*sizeX] - w_n[x];
            w[x + y*sizeX] = w_n[x];
        }

        // swap history
        double *tmpSwap = w_n3; // remember where we allocated that piece of memory, so we can reuse it in the next iteration
        w_n3 = w_n2;
        w_n2 = w_n1;
        w_n1 = w_n;
        w_n = tmpSwap;

    }
}


// ---- filter backward direction y-axis

void AlgAnisotropicGaussBase::filterYoungBackwardY(double* lastLineOfInput, double* w, double* output, filterLines outputLines, coefficients coeffs)
{
    double *tmpSwap;
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};

    double* output_n = outputLines.line_n;
    double* output_n1 = outputLines.line_n1;
    double* output_n2 = outputLines.line_n2;
    double* output_n3 = outputLines.line_n3;

    // prepare values for first iteration (incl Triggs)
    filterYoungBackwardInitialValuesY(lastLineOfInput, output, outputLines, coeffs);

    // iteration
    double sqrtb4 = sqrt(b[4]);
    for (int y = sizeY - 2; y >= 0; y--)
    {
        output_n[sizeX - 1] = b[1]*output_n1[sizeX - 1] + b[2]*output_n2[sizeX - 1] + b[3]*output_n3[sizeX - 1];
        output_n[sizeX - 1] = sqrtb4*w[sizeX - 1 + y*sizeX] - output_n[sizeX - 1];
        output[sizeX - 1 + y*sizeX] = output_n[sizeX - 1];

        for (int x = sizeX - 2; x >= 0; x--)
        {
            // step
            output_n[x] = b[1]*output_n1[x] + b[2]*output_n2[x] + b[3]*output_n3[x];
            output_n[x] = sqrtb4*w[x + y*sizeX] - output_n[x];
            output[x + y*sizeX] = output_n[x];
        }

        // swap history
        tmpSwap = output_n3; // remember where we allocated that piece of memory, so we can reuse it in the next iteration
        output_n3 = output_n2;
        output_n2 = output_n1;
        output_n1 = output_n;
        output_n = tmpSwap;
    }

}


void AlgAnisotropicGaussBase::filterYoungBackwardInitialValuesY(double *lastLineOfInput, double* output, filterLines outputLines, coefficients coeffs)
{
    int sizeX = m_sizeX;
    int sizeY = m_sizeY;
    double m11 = coeffs.M[0], m12 = coeffs.M[1], m13 = coeffs.M[2];
    double m21 = coeffs.M[3], m22 = coeffs.M[4], m23 = coeffs.M[5];
    double m31 = coeffs.M[6], m32 = coeffs.M[7], m33 = coeffs.M[8];
    double b[5] = {coeffs.b[0], coeffs.b[1], coeffs.b[2], coeffs.b[3], coeffs.b[4]};

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

    double sqrtb4 = sqrt(b[4]);


    for (int x = sizeX - 1; x >= 0; x--)
    {
        u_plus = lastLineOfInput[x]*sqrtb4/(1.0 - b1 - b2 - b3);
        v_plus = u_plus/(1.0 - b1 - b2 - b3);

        u_t = output_n1[x] - u_plus;
        u_t1 = output_n2[x] - u_plus;
        u_t2 = output_n3[x] - u_plus;

        output_n1[x] = (m11*u_t + m12*u_t1 + m13*u_t2 + v_plus)*sqrtb4;
        output[(sizeY - 1)*sizeX + x] = output_n1[x];
        output_n2[x] = (m21*u_t + m22*u_t1 + m23*u_t2 + v_plus)*sqrtb4;
        output_n3[x] = (m31*u_t + m32*u_t1 + m33*u_t2 + v_plus)*sqrtb4;

    }


}



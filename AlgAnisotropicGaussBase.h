#ifndef ANIGAUSSHYBBASE_H
#define ANIGAUSSHYBBASE_H
#include <System/SizedIntegers.h>
#include "structs.h"

/**
 * @class AlgAnisotropicGaussBase.h
 * @brief anisotropic Gaussian filter (2D)
 *
 * The image is recursively filtered with an anisotropic Gaussian filter in 2D
 * which is determined by the standard deviations sigma_u, sigma_v and
 * a rotation angle theta.
 *
 * This is the base class, so other classes can inherit from it. The other classes
 * are different w.r.t. the filter direction and the interpolation, which only happens
 * in t-direction, i.e. the t-line.
 *
 *
 * @bib
 *      - J.-M. Geusebroek, A. Smeulders, and J. Weijer.
 *        "Fast anisotropic gauss filtering." Volume 2350, 09 2003.
 *      - B. Triggs and M. Sdika. "Boundary conditions for young-van vliet recursive
 *        filtering." Signal Processing, IEEE Transactions on, 54:2365 – 2367,
 *        07 2006.
 *      - I. Young, L. Van Vliet, and M. van Ginkel. "Recursive gabor filtering."
 *        Signal Processing, IEEE Transactions on, 50:2798 – 2805, 12 2002.
 */
class AlgAnisotropicGaussBase
{
public:
    /**
     * @brief Constructor.
     * @param input[in] - pointer to the input image
     * @param sizeX[in] - size of images wrt x-dimension
     * @param sizeY[in] - size of images wrt y-dimension
     */
    AlgAnisotropicGaussBase(double *input, int sizeX, int sizeY);
    ~AlgAnisotropicGaussBase();

    virtual void runFilter(double *output, double sigma_u, double sigma_v, double theta){}

protected:
    void computeDeviationsAndIntercept(double sigma_u, double sigma_v, double theta, bool flipXY);
    double calculateQ(double sigma);
    void calculateBCoefficients(double q, coefficients* coeffs);
    void calculateM(coefficients* coeffs);

    // filter y-axis
    void filterXAxis(double *input, double *w, double *output, double sigma);
    void filterYoungForwardX(double *input,  double* w, coefficients coeffs);
    void filterYoungBackwardX(double* input, double *lastColumnOfInput, double* w, double* output, coefficients coeffs);
    void filterYoungBackwardInitialValuesX(double* input, double *lastColumnOfInput, double* w, int offset, coefficients coeffs);


    // filter x-axis
    void filterYAxis(double *input, double *w, double *output, double sigma);
    void filterYoungForwardY(double *input, double* w, filterLines wLines, coefficients coeffs);
    void filterYoungBackwardY(double* lastLineOfInput, double* w, double* output, filterLines outputLines, coefficients coeffs);
    void filterYoungBackwardInitialValuesY(double* lastLineOfInput, double* output, filterLines outputLines, coefficients coeffs);

       //variables
    int m_sizeX;
    int m_sizeY;

    double m_sigma_x;
    double m_sigma_phi;
    double m_tan_phi;

    double* m_input;

    double* m_initialValuesBackwards;

};

#endif // ANIGAUSSHYBBASE_H

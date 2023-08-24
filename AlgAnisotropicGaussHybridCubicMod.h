#ifndef ANIGAUSSCUBICS_H
#define ANIGAUSSCUBICS_H
#include <System/SizedIntegers.h>
#include "structs.h"
#include "AlgAnisotropicGaussBase.h"


/**
 * @class AlgAnisotropicGaussHybridCubicMod.h
 * @brief anisotropic Gaussian filter (2D) hybrid with cubic interpolation +
 * major-axis modification
 *
 * The image is recursively filtered with an anisotropic Gaussian filter in 2D
 * which is determined by the standard deviations sigma_u, sigma_v and
 * a rotation angle theta.
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
class AlgAnisotropicGaussHybridCubicMod: protected AlgAnisotropicGaussBase
{
public:
    /**
     * @brief Constructor.
     * @param input[in] - pointer to the input image
     * @param sizeX[in] - size of images wrt x-dimension
     * @param sizeY[in] - size of images wrt y-dimension
     */
    AlgAnisotropicGaussHybridCubicMod(double *input, int sizeX, int sizeY);
    virtual ~AlgAnisotropicGaussHybridCubicMod();

    /**
     * @brief The anisotropic Gaussian filter is run.
     * @param output[out] - pointer to the output image
     * @param sigma_u[in] - standard deviation of the Gaussian for the major axis
     * @param sigma_v[in] - standard deviation of the Gaussian for the minor axis
     * @param theta[in] - rotation angle
     */
    void runFilter(double *output, double sigma_u, double sigma_v, double theta);

private:

    void calculateY2(int x, double* y, double* y2, double *u);
    void calculateX2(int x, double* y, double* y2, double *u);

    //filter t-line
    void filterTLine(double *input, double *w, double *output);
    void filterYoungForwardT(double *input, double* w, filterLines wLines, coefficients coeffs);
    void filterYoungBackwardT(double* input, double* w, double* output, filterLines outputLines, coefficients coeffs);
    void filterYoungBackwardInitialValuesT(double* lastLineOfInput, double* output, filterLines outputLines, coefficients coeffs);

    //filter t-line along x
    void filterTxLine(double *input, double *w, double *output);
    void filterYoungForwardTx(double *input, double* w, filterLines wLines, coefficients coeffs);
    void filterYoungBackwardTx(double* lastLineOfInput, double* w, double* output, filterLines outputLines, coefficients coeffs);
    void filterYoungBackwardInitialValuesTx(double* lastLineOfInput, double* output, filterLines outputLines, coefficients coeffs);

};

#endif // ANIGAUSSCUBICS_H

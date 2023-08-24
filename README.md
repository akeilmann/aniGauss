# aniGauss
AniGauss is a C++ implementation for anisotropic Gaussian filters in 2D.
This repository provides the hybrid algorithm with both linear and cubic implementation, as well as with and without modification.

You can contact us via [keilmann@rptu.de](mailto:keilmann@rptu.de).
    
## Publication
For details, please see the paper

[A. Keilmann, M. Godehardt, A. Moghiseh, C. Redenbach, and K. Schladitz,
     ‘Improved Anisotropic Gaussian Filters’, arXiv [eess.IV]. 2023.](https://arxiv.org/pdf/2303.13278.pdf).

If you use our code and write a paper, please cite us:

```
@misc{keilmann2023improved,
      title={Improved Anisotropic Gaussian Filters}, 
      author={Alex Keilmann and Michael Godehardt and Ali Moghiseh and Claudia Redenbach and Katja Schladitz},
      year={2023},
      eprint={2303.13278},
      archivePrefix={arXiv},
      primaryClass={eess.IV}
}
```

Further, this work is based on the following papers:

- J.-M. Geusebroek, A. Smeulders, and J. Weijer. "Fast anisotropic gauss filtering." Volume 2350, 09 2003.

- B.&nbsp; Triggs and M. Sdika. "Boundary conditions for young-van vliet recursive filtering." Signal Processing, IEEE Transactions on, 54:2365 – 2367, 07 2006.

- I.&nbsp; Young, L. Van Vliet, and M. van Ginkel. "Recursive gabor filtering." Signal Processing, IEEE Transactions on, 50:2798 – 2805, 12 2002.

## Code Example

This example applies an anisotropic Gaussian filter with sigma_u = 2.0, sigma_v = 25.0 and phi = 60° to an image (unit impulse).
It uses linear interpolation and the modification proposed in the paper (see above).

```
int nSizeX = 256, nSizeY = 256;
double sigma_u = 2.0, sigma_v = 25.0;
double phi = 60;

double *pImgInDouble = new double[nSizeX * nSizeY];
memset(pImgInDouble, 0, nSizeX * nSizeY * sizeof(*pImgInDouble);
pImgInDouble[nSizeX/2.0*nSizeY + nSizeX/2.0] = 1.0;

double *pImgTmp = new double[nSizeX * nSizeY];

AlgAnisotropicGaussHybridLinearMod* algHybridLinearMod = 
    new AlgAnisotropicGaussHybridLinearMod(pImgInDouble, nSizeX, nSizeY);
algHybridLinearMod->runFilter(pImgTmp, sigma_u, sigma_v, phi);
```

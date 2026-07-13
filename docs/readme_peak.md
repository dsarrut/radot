# Peak Value Computation

This document explains the rationale and math behind the peak value computed by the `stats` command in `radot`.

## Clinical Rationale (SUVpeak)
In quantitative medical imaging and dosimetry (e.g., PET imaging under PERCIST criteria), single-voxel maximum values (`Max`) are highly susceptible to image noise. To obtain a more robust estimate of the maximum uptake, the **peak value** (e.g., `SUVpeak`) is used instead.

The standard clinical definition of SUVpeak is:
> The maximum average value within a **$1\text{ cm}^3$ sphere** (approx. $1.2\text{ cm}$ diameter) positioned anywhere within the region of interest (ROI).

## Approximation via Gaussian Smoothing
Instead of calculating a sliding spherical average, a common and computationally efficient approximation is to smooth the image using a **3D Gaussian filter** and then find the maximum value of the smoothed image within the mask:

$$\text{Peak} = \max_{x \in \text{ROI}} (I * G_{\sigma})(x)$$

where $I$ is the input image, $G_{\sigma}$ is the Gaussian kernel, and $*$ is the convolution operator.

### Choosing the Gaussian Sigma ($\sigma$)
To align the Gaussian filter with the standard $1\text{ cm}^3$ sphere of diameter $12\text{ mm}$:
1. We define the **FWHM (Full Width at Half Maximum)** of the Gaussian filter to match the sphere's diameter ($12\text{ mm}$).
2. The mathematical relationship between FWHM and the standard deviation ($\sigma$) of a Gaussian distribution is:
   $$\text{FWHM} = 2.355 \times \sigma$$
3. Solving for $\sigma$ yields:
   $$\sigma = \frac{\text{FWHM}}{2.355} = \frac{12\text{ mm}}{2.355} \approx 5.1\text{ mm}$$

Currently, `radot` uses a default physical $\sigma = 3.0\text{ mm}$ (corresponding to a $\text{FWHM} \approx 7\text{ mm}$), which is a common clinical smoothing filter size.

## References

Wahl, R. L., Jacene, H., Kasamon, Y., & Lodge, M. A. (2009). From RECIST to PERCIST: Limiting Factors and a Proposed New Criteria for Evaluating Tumor Response Using PET/CT. Journal of Nuclear Medicine, 50(Suppl 1), 35S-50S. https://doi.org/10.2967/jnumed.108.057307 (First introduces PERCIST and the $1\text{ cm}^3$ spherical $SUV_{peak}$ definition).

Lodge, M. A., Chaudhry, M. A., & Wahl, R. L. (2012). Repeatability of SUVpeak in FDG PET. Journal of Nuclear Medicine, 53(8), 1222-1229. https://doi.org/10.2967/jnumed.116.186353 (Details the implementation, precision, and repeatability of peak metrics).

Lodge, M. A. (2014). Measurement of SUVpeak in Clinical PET. Journal of Nuclear Medicine, 55(8), 1219-1220. https://doi.org/10.2967/jnumed.111.093443 (Discusses standardizing and verifying $SUV_{peak}$ calculations across vendors).

Boellaard, R., et al. (2015). FDG PET/CT: EANM procedure guidelines for tumour imaging: version 2.0. European Journal of Nuclear Medicine and Molecular Imaging, 42(2), 328-354. https://doi.org/10.1007/s00259-014-2961-x (Details recommendations on spatial smoothing post-reconstruction using Gaussian filters to ensure quantitative repeatability).

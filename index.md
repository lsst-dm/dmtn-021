#Implementation of Image Difference Decorrelation

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

## Abstract

Herein, we describe ...

## Introduction

The standard method for PSF matching for image subtraction in the LSST software stack is the method of [Alard & Lupton (1998)](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) (hereafter *A&L*). This algorithm learns a convolution kernel which, when convolved with the template image, matches the PSF of the template with that of the science image. Due to its use of linear basis functions to model the matching kernel, the method can cleanly incorporate spatially-varying PSFs (i.e., via a spatially-varying matching kernel), as well as a spatially-varying differential background. The algorithm has the advantage that it does not require measurement of the images' PSFs. Instead it only needs to model the differential (potentially spatially-varying) matching kernel in order to obtain an accurate subtraction.

The [Alard & Lupton (1998)](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) method produces an optimal difference image in the case of a noise-less template. However, when the template is noisy (e.g., when the template is comprised of a small number of co-adds), then its convolution with the matching kernel leads to significant covariance of neighboring pixels within the subtracted image, which will affect detection and measurement if not accounted for ([Slater, et al. 2016](http://dmtn-006.lsst.io)). False detections in this case can be reduced by tracking the covariance matrix, or more *ad-hoc* (as is the current implementation) increasing the detection threshold. 

While LSST will, over its ten-year span, collect dozens of observations per field and passband, at the onset of the survey, this number will be small enough that this issue of noisy templates will be important. Moreover, if we inted to bin templates by airmass to account for differential chromatic refraction (DCR), then the total number of coadds contributing to each template will by necessity be smaller. Finally, depending upon the flavor of coadd ([Bosch, 2016](http://dmtn-015.lsst.io)) used to construct the template, template noise and the resulting covariances in the image difference will be more or less of an issue as the survey progresses.

## Proposal

An algorithm developed by [Kaiser, 2001](Addition of Images with Varying Seeing. PSDC-002-011-xx) and later rediscovered by [Zackay, et al (2015)](https://arxiv.org/abs/1512.06879) showed that the noise in a PSF-matched coadd image can be decorrelated via noise whitening (i.e. flattening the noise spectrum). The same principle may also be applied to image differencing ([Zackay, et al. (2016)](https://arxiv.org/abs/1601.02655)). In the case of [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) - based PSF matching, this results in an image difference in Fourier space $\widehat{D}(k)$ (where $\widehat{x}(k)$ denotes the Fourier transform of $D$): 

###### Equation 1.
\$\$
	\\widehat{D}(k) = [\\widehat{I}_1(k) - \\widehat{\\kappa}(k) \\widehat{I}_2(k)] \\sqrt{ \\frac{ \\sigma_1^2 + \\sigma_2^2}{ \\sigma_1^2 + \\widehat{\\kappa}^2(k) \\sigma_2^2}}
\$\$

Here, $I_1$ and $I_2$ are the two images being subtracted (typically $I_2$ is the template image, which is convolved with the PSF matching kernel $\kappa$). $\sigma_1^2$ and $\sigma_2^2$ are the variances of the two respective images. The term in the square-root of is a *post-subtraction convolution kernel* $\widehat{\phi}(k)$, which, when convolved with the image difference, has the effect of decorrelating the noise in the image difference. Thus, we may perform PSF matching to estimate $\kappa$ by standard methods (e.g., [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) and related methods) and then correct for the noise in the template. This maintains the advantages described previously: the PSFs of $I_1$ and $I_2$ do not need to be measured, and spatial variations in PSFs may be readily accounted for (although see below). The decorrelation can be relatively inexpensive, as it requires (at least) one *FFT* of $\kappa$ and *iFFT* of $\widehat{\phi}(k)$ (which are both small, of the order 1,000 pixels), followed by one convolution.

## Implementation details

Since the current implementation of [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) is performed in image space, we chose to implement the image decorrelation in image space as well. The image differencing is performed as usual to estimate $\kappa$ and compute the uncorrected image difference, $I_1 - (\kappa \otimes I_2)$. The *post-subtraction convolution kernel* $\widehat{\phi}(k)$ is then computed in frequency space from $\widehat{\kappa}(k)$, $\sigma_1$, and $\sigma_2$, and is then inverse Fourier-transformed to a kernel $\phi$ in real space. The image difference is then convolved with $\phi$ to obtain the decorrelated image difference, $D(x) = \phi \otimes \big[ I_1 - (\kappa \otimes I_2) \big]$. 

## Results

## Conclusions and future work

### Accounting for spatial variations in PSF matching kernel and noise



## References

## Appendix
### Appendix A. Implementation of basic Zackay et al. (2016) algorithm.
### Appendix B. Something else.

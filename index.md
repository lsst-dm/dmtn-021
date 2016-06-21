#Implementation of Image Difference Decorrelation

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

# Abstract

Herein, we describe ...

# Introduction

The standard method for PSF matching for image subtraction in the LSST software stack is the method of [Alard & Lupton (1998)](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) (hereafter *A&L*). This algorithm learns a convolution kernel which, when convolved with the template image, matches the PSF of the template with that of the science image. Due to its use of linear basis functions to model the matching kernel, the method can cleanly incorporate spatially-varying PSFs (i.e., via a spatially-varying matching kernel), as well as a spatially-varying differential background. The algorithm has the advantage that it does not require measurement of the images' PSFs. Instead it only needs to model the differential (potentially spatially-varying) matching kernel in order to obtain an accurate subtraction.

The [Alard & Lupton (1998)](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) method produces an optimal difference image in the case of a noise-less template. However, when the template is noisy (e.g., when the template is comprised of a small number of co-adds), then its convolution with the matching kernel leads to significant covariance of neighboring pixels within the subtracted image, which will affect detection and measurement if not accounted for ([Slater, et al. 2016](http://dmtn-006.lsst.io)). False detections in this case can be reduced by tracking the covariance matrix, or more *ad-hoc* (as is the current implementation) increasing the detection threshold. 

While LSST will, over its ten-year span, collect dozens of observations per field and passband, at the onset of the survey, this number will be small enough that this issue of noisy templates will be important. Moreover, if we inted to bin templates by airmass to account for differential chromatic refraction (DCR), then the total number of coadds contributing to each template will by necessity be smaller. Finally, depending upon the flavor of coadd ([Bosch, 2016](http://dmtn-015.lsst.io)) used to construct the template, template noise and the resulting covariances in the image difference will be more or less of an issue as the survey progresses.

# Proposal

An algorithm developed by [Kaiser, 2001](Addition of Images with Varying Seeing. PSDC-002-011-xx) and later rediscovered by [Zackay, et al (2015)](https://arxiv.org/abs/1512.06879) showed that the noise in a PSF-matched coadd image can be decorrelated via noise whitening (i.e. flattening the noise spectrum). The same principle may also be applied to image differencing ([Zackay, et al. (2016)](https://arxiv.org/abs/1601.02655)). In the case of [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) - based PSF matching, this results in an image difference in Fourier space $\widehat{D}(k)$ (where $\widehat{x}(k)$ denotes the Fourier transform of $D$): 

$$
\widehat{D}(k) = \big[ \widehat{I}_1(k) - \widehat{\kappa}(k) \widehat{I}_2(k) \big] \sqrt{ \frac{ \sigma_1^2 + \sigma_2^2}{ \sigma_1^2 + \widehat{\kappa}^2(k) \sigma_2^2}}
$$

######Equation 1.

Here, $I_1$ and $I_2$ are the two images being subtracted (typically $I_2$ is the template image, which is convolved with the PSF matching kernel $\kappa$). $\sigma_1^2$ and $\sigma_2^2$ are the variances of the two respective images. The term in the square-root of is a *post-subtraction convolution kernel* $\widehat{\phi}(k)$, which, when convolved with the image difference, has the effect of decorrelating the noise in the image difference. Thus, we may perform PSF matching to estimate $\kappa$ by standard methods (e.g., [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) and related methods) and then correct for the noise in the template. This maintains the advantages described previously: the PSFs of $I_1$ and $I_2$ do not need to be measured, and spatial variations in PSFs may be readily accounted for (although see below). The decorrelation can be relatively inexpensive, as it requires (at least) one *FFT* of $\kappa$ and *iFFT* of $\widehat{\phi}(k)$ (which are both small, of the order 1,000 pixels), followed by one convolution.

# Implementation details

Since the current implementation of [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) is performed in image space, we chose to implement the image decorrelation in image space as well. The image differencing is performed as usual to estimate $\kappa$ and compute the uncorrected image difference, $I_1 - (\kappa \otimes I_2)$. The *post-subtraction convolution kernel* $\widehat{\phi}(k)$ is then computed in frequency space from $\widehat{\kappa}(k)$, $\sigma_1$, and $\sigma_2$, and is then inverse Fourier-transformed to a kernel $\phi$ in real space. The image difference is then convolved with $\phi$ to obtain the decorrelated image difference, $D(x) = \phi \otimes \big[ I_1 - (\kappa \otimes I_2) \big]$. 

# Results

We have developed a simple reference implementation of [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A), and applied it to simulated images with point-sources with a variety of signal-to-noise, and different Gaussian PSFs and image variances. We included the capability to simulate spatial PSF variation, including spatially-varying astrometric offsets (which can be incorporated into the [A&L](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) PSF matching kernel). An example input template and science image, as well as PSF-matched template and resulting *diffim* is shown in [Figure 1](#figure-1-image-differencing).

In [Figure 2](#figure-2-kernels), we show the PSF matching kernel ($\kappa$) that was estimated for the images shown in [Figure 1](#figure-1-image-differencing), and the resulting decorrelation kernel, $\phi$. We note that $\phi$ largely has the structure of a delta function, with a small region of negative signal, thus its capability, when convolved with the difference image, to act as an effective "sharpening" kernel.

![](_static/img0.png)

<a name="figure-1-image-differencing"/></a>

###### *Figure 1. Image differencing.*

*From left to right, sample (simulated) template image, PSF-matched template, science image, and difference image. In this simulated example, the source near the center was set to increase in flux by 2% between the science and template "exposures."*

![Matching kernel](_static/img1.png)
![Correction kernel](_static/img2.png)

<a name="figure-2-kernels"/></a>

###### *Figure 2. Kernels.*

*Sample PSF matching kernel $\kappa$ (left) and resulting decorrelation kernel, $\phi$ for the images shown in [Figure 1](#figure-1-image-differencing).*

When we convolve $\phi$ ([Figure 2](#figure-2-kernels), right panel) with the raw image difference ([Figure 1](#figure-1-image-differencing), right-most panel), we obtain the decorrelated image, shown in the left-most panel of [Figure 3](#figure-3-decorrelated-diffim). While the noise visually appears to be greater in the decorrelated image, a closer look at the statistics reveals that this is indeed the case ([Figure 4](#figure-4-decorrelated-image-statistics) and [Figure 5](figure-5)). [Figure 4](#figure-4-decorrelated-image-statistics) shows that the variance of the decorrelated image has increased. Indeed, the measured variances reveal that the variance of the uncorrected image difference was lower than expected, while the decorrelation has increased the variance to the expected level:

```python
%In [1]:
print sig1, sig2  # Input std. deviation of template and science images
print 'Corrected:', np.mean(diffim2), np.std(diffim2)
print 'Original: ', np.mean(diffim1), np.std(diffim1)
print 'Expected: ', np.sqrt(sig1**2 + sig2**2)
%Out [1]:
0.2 0.2
Corrected: 10.0042330181 0.293237231242
Original:  9.99913482654 0.211891941431
Expected:  0.282842712475
```

In addition, we see ([Figure 5](#figure-5-covariance-matrices)) that the covariances between neighboring pixels in the image difference has been significantly decreased following convolution with the decorrelation kernel. The covariance matrix has been significantly diagonalized:

```python
%In [2]:
print np.nansum(cov2)/np.sum(np.diag(cov2))  # cov2 is the covar. matrix of the corrected image.
print np.nansum(cov1)/np.sum(np.diag(cov1))  # cov1 is the covar. matrix of the uncorrected image.
%Out [2]:
0.300482626371
0.793176605206
```

![Decorrelated diffim](_static/img3.png)

<a name="figure-3-decorrelated-diffim"/></a>

###### *Figure 3. Decorrelated diffim.*

*On the left is the decorrelated image difference. Original image difference is shown here for comparison, in the right-most panel, with the same intensity scale, as well as in* [Figure 1](#figure-1-image-differencing).

![Decorrelated image statistics](_static/img4.png)

<a name="figure-4-decorrelated-image-statistics"/></a>

###### *Figure 4. Decorrelated image statistics.*

*Histogram of sigma-clipped pixels in the original image difference (blue; 'orig') and the decorrelated image difference (red; 'corr') in* [Figure 3](#figure-3-decorrelated-diffim).

![Covariance matrix 1](_static/img5.png)
![Covariance matrix 2](_static/img6.png)

<a name="figure-5-covariance-matrices"/></a>

###### *Figure 5. Covariance matrices.*

*Covariance between neighboring pixels in the original, uncorrected image difference (left) and the decorrelated image difference (right) in* [Figure 3](#figure-3-decorrelated-diffim).

# Conclusions and future work

### Accounting for spatial variations in noise and matching kernel

# References

# Appendix

### Appendix A. Implementation of basic Zackay et al. (2016) algorithm.

### Appendix B. Notebooks and code

All figures in this document and related code are from notebooks in [the diffimTests github repository](https://github.com/lsst-dm/diffimTests), in particular, [this](https://github.com/djreiss/diffimTests/blob/master/14.%20Test%20Lupton(ZOGY)%20post%20convolution%20kernel%20on%20simulated%20(noisy)%202-D%20data%20with%20a%20variable%20source-updated.ipynb) and [this](https://github.com/djreiss/diffimTests/blob/master/13.%20compare%20L(ZOGY)%20and%20ZOGY%20diffims%20and%20PSFs.ipynb) one.

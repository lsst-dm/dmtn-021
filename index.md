#Implementation of Image Difference Decorrelation

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

## Abstract

Herein, we describe ...

## Introduction

The standard method for PSF matching for image subtraction in the LSST software stack is the method of [Alard & Lupton (1998)](http://adsabs.harvard.edu/abs/1998ApJ...503..325A). This algorithm learns a convolution kernel which, when convolved with the template image, matches the PSF of the template with that of the science image. Do to its use of linear basis functions to model the matching kernel, the method can cleanly incorporate spatially-varying PSFs (i.e., via a spatially-varying matching kernel), as well as a spatially-varying differential background. The algorithm has the advantage that it does not require measurement of the images' PSFs. Instead it only needs to model the differential (potentially spatially-varying) matching kernel in order to obtain an accurate subtraction.

The [Alard & Lupton (1998)](http://adsabs.harvard.edu/abs/1998ApJ...503..325A) method produces an optimal difference image in the case of a noise-less template. However, when the template is noisy (e.g., when the template is comprised of a small number of co-adds), then its convolution with the matching kernel leads to significant covariance of neighboring pixels within the subtracted image, which will affect detection if not accounted for ([Slater, et al. 2016](http://dmtn-006.lsst.io)). False detections in this case can be reduced by tracking the covariance, or more simply (as is the current implementation) increasing the detection threshold. 

While LSST will, over its ten-year span, collect dozens of observations per field and passband, at the onset of the survey, this number will be small enough that this issue of noisy templates will be important. Moreover, if we inted to bin templates by airmass to account for differential chromatic refraction (DCR), then the total number of coadds contributing to each template will by necessity be smaller. Finally, depending upon the flavor of coadd ([Bosch, 2016](http://dmtn-015.lsst.io)) used to construct the template, template noise and the resulting covariances in the image difference will be more or less of an issue as the survey progresses.

## Proposal

\$\$
	D(k) = [I_1(k) - \\kappa(k) I_2(k)] \\sqrt{ \\frac{ \\sigma_1^2 + \\sigma_2^2}{ \\sigma_1^2 + \\kappa^2(k) \\sigma_2^2}}
\$\$

###### Equation 1.

## Implementation

## Conclusions and future work

## References

## Appendix
### Appendix A. Implementation of basic Zahav et al. (2016) algorithm.
### Appendix B. Something else.

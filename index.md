#Implementation of Image Difference Decorrelation

## Abstract



## Introduction

The standard method for PSF matching to produce image differences in the LSST software stack is the method of Alard & Lupton (1998). This algorithm learns a convolution kernel which, when the template image is convolved with it, matches the PSF of the template with the science image. Do to its use of linear basis functions to model the matching kernel, the method can cleanly model a spatially-varying matching kernel, as well as a spatially-varying differential background. The algorithm has the advantage that it does not require accurate measurements of the images' PSFs. Instead it only needs to model the differential (potentially spatially-varying) matching kernel in order to obtain an accurate subtraction.

A disadvantage of this method is that if the template contains significant noise (which will be the case when the template is comprised of a small number of co-adds) is that the convolution of the template by the matching kernel leads to significant covariance of neighboring pixels within the subtracted image, which will affect detection (http://dmtn-006.lsst.io).

## Proposal

## Implementation

## Conclusions and future work

## References

## Appendix
### Appendix A. Implementation of basic Zahav et al. (2016) algorithm.
### Appendix B. Something else.

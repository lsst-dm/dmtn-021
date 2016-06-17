Implementation of Image Difference Decorrelation
================================================

Abstract
--------

Herein, we describe ...

Introduction
------------

The standard method for PSF matching for image subtraction in the LSST
software stack is the method of Alard & Lupton (1998). This algorithm
learns a convolution kernel which, when convolved with the template
image, matches the PSF of the template with that of the science image.
Do to its use of linear basis functions to model the matching kernel,
the method can cleanly incorporate spatially-varying PSFs (i.e., via a
spatially-varying matching kernel), as well as a spatially-varying
differential background. The algorithm has the advantage that it does
not require measurement of the images' PSFs. Instead it only needs to
model the differential (potentially spatially-varying) matching kernel
in order to obtain an accurate subtraction.

The Alard & Lupton (1998) method produces an optimal difference image in
the case of a noise-less template. However, when the template is noisy
(e.g., when the template is comprised of a small number of co-adds),
then its convolution with the matching kernel leads to significant
covariance of neighboring pixels within the subtracted image, which will
affect detection if not accounted for (http://dmtn-006.lsst.io).

Proposal
--------

Implementation
--------------

Conclusions and future work
---------------------------

References
----------

Appendix
--------

Appendix A. Implementation of basic Zahav et al. (2016) algorithm.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Appendix B. Something else.
~~~~~~~~~~~~~~~~~~~~~~~~~~~

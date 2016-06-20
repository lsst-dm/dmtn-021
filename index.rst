:tocdepth: 1
Implementation of Image Difference Decorrelation
================================================

.. raw:: html

   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

Abstract
========

Herein, we describe ...

Introduction
============

The standard method for PSF matching for image subtraction in the LSST
software stack is the method of `Alard & Lupton
(1998) <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`__ (hereafter
*A&L*). This algorithm learns a convolution kernel which, when convolved
with the template image, matches the PSF of the template with that of
the science image. Due to its use of linear basis functions to model the
matching kernel, the method can cleanly incorporate spatially-varying
PSFs (i.e., via a spatially-varying matching kernel), as well as a
spatially-varying differential background. The algorithm has the
advantage that it does not require measurement of the images' PSFs.
Instead it only needs to model the differential (potentially
spatially-varying) matching kernel in order to obtain an accurate
subtraction.

The `Alard & Lupton
(1998) <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`__ method
produces an optimal difference image in the case of a noise-less
template. However, when the template is noisy (e.g., when the template
is comprised of a small number of co-adds), then its convolution with
the matching kernel leads to significant covariance of neighboring
pixels within the subtracted image, which will affect detection and
measurement if not accounted for (`Slater, et al.
2016 <http://dmtn-006.lsst.io>`__). False detections in this case can be
reduced by tracking the covariance matrix, or more *ad-hoc* (as is the
current implementation) increasing the detection threshold.

While LSST will, over its ten-year span, collect dozens of observations
per field and passband, at the onset of the survey, this number will be
small enough that this issue of noisy templates will be important.
Moreover, if we inted to bin templates by airmass to account for
differential chromatic refraction (DCR), then the total number of coadds
contributing to each template will by necessity be smaller. Finally,
depending upon the flavor of coadd (`Bosch,
2016 <http://dmtn-015.lsst.io>`__) used to construct the template,
template noise and the resulting covariances in the image difference
will be more or less of an issue as the survey progresses.

Proposal
========

An algorithm developed by `Kaiser,
2001 <Addition%20of%20Images%20with%20Varying%20Seeing.%20PSDC-002-011-xx>`__
and later rediscovered by `Zackay, et al
(2015) <https://arxiv.org/abs/1512.06879>`__ showed that the noise in a
PSF-matched coadd image can be decorrelated via noise whitening (i.e.
flattening the noise spectrum). The same principle may also be applied
to image differencing (`Zackay, et al.
(2016) <https://arxiv.org/abs/1601.02655>`__). In the case of
`A&L <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`__ - based PSF
matching, this results in an image difference in Fourier space
:math:`\widehat{D}(k)` (where :math:`\widehat{x}(k)` denotes the Fourier
transform of :math:`D`):

.. raw:: html

   <table>

.. math::


   \widehat{D}(k) = \big[ \widehat{I}_1(k) - \widehat{\kappa}(k) \widehat{I}_2(k) \big] \sqrt{ \frac{ \sigma_1^2 + \sigma_2^2}{ \sigma_1^2 + \widehat{\kappa}^2(k) \sigma_2^2}}

Equation 1.

.. raw:: html

   </table>

Here, :math:`I_1` and :math:`I_2` are the two images being subtracted
(typically :math:`I_2` is the template image, which is convolved with
the PSF matching kernel :math:`\kappa`). :math:`\sigma_1^2` and
:math:`\sigma_2^2` are the variances of the two respective images. The
term in the square-root of is a *post-subtraction convolution kernel*
:math:`\widehat{\phi}(k)`, which, when convolved with the image
difference, has the effect of decorrelating the noise in the image
difference. Thus, we may perform PSF matching to estimate :math:`\kappa`
by standard methods (e.g.,
`A&L <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`__ and related
methods) and then correct for the noise in the template. This maintains
the advantages described previously: the PSFs of :math:`I_1` and
:math:`I_2` do not need to be measured, and spatial variations in PSFs
may be readily accounted for (although see below). The decorrelation can
be relatively inexpensive, as it requires (at least) one *FFT* of
:math:`\kappa` and *iFFT* of :math:`\widehat{\phi}(k)` (which are both
small, of the order 1,000 pixels), followed by one convolution.

Implementation details
======================

Since the current implementation of
`A&L <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`__ is performed
in image space, we chose to implement the image decorrelation in image
space as well. The image differencing is performed as usual to estimate
:math:`\kappa` and compute the uncorrected image difference,
:math:`I_1 - (\kappa \otimes I_2)`. The *post-subtraction convolution
kernel* :math:`\widehat{\phi}(k)` is then computed in frequency space
from :math:`\widehat{\kappa}(k)`, :math:`\sigma_1`, and
:math:`\sigma_2`, and is then inverse Fourier-transformed to a kernel
:math:`\phi` in real space. The image difference is then convolved with
:math:`\phi` to obtain the decorrelated image difference,
:math:`D(x) = \phi \otimes \big[ I_1 - (\kappa \otimes I_2) \big]`.

Results
=======

We have developed a simple reference implementation of
`A&L <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`__, and applied
it to simulated images with point-sources with a variety of
signal-to-noise, and different Gaussian PSFs and image variances. We
included the capability to simulate spatial PSF variation, including
spatially-varying astrometric offsets (which can be incorporated into
the `A&L <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`__ PSF
matching kernel). An example input template and science image, as well
as PSF-matched template and resulting *diffim* is shown in `Figure
1 <#figure-1-image-differencing>`__.

In `Figure 2 <#figure-2-kernels>`__, we show the PSF matching kernel
(:math:`\kappa`) that was estimated for the images shown in `Figure
1 <#figure-1-image-differencing>`__, and the resulting decorrelation
kernel, :math:`\phi`. We note that :math:`\phi` largely has the
structure of a delta function, with a small region of negative signal,
thus its capability, when convolved with the difference image, to act as
an effective "sharpening" kernel.

.. figure:: _static/img0.png
   :alt: 

Figure 1. Image differencing.
-----------------------------

From left to right, sample (simulated) template image, PSF-matched
template, science image, and difference image. In this simulated
example, the source near the center was set to increase in flux by 2%
between the science and template "exposures."

|Matching kernel| |Correction kernel|

Figure 2. Kernels.
------------------

Sample PSF matching kernel :math:`\kappa` (left) and resulting
decorrelation kernel, :math:`\phi` for the images shown in `Figure
1 <#figure-1-image-differencing>`__.

Conclusions and future work
===========================

Accounting for spatial variations in noise and matching kernel
--------------------------------------------------------------

References
==========

Appendix
========

Appendix A. Implementation of basic Zackay et al. (2016) algorithm.
-------------------------------------------------------------------

Appendix B. Something else.
---------------------------

.. |Matching kernel| image:: _static/img1.png
.. |Correction kernel| image:: _static/img2.png

This is a MatLab/C++ library of our ICCV 2013 paper [1].

![ScreenShot](screenshot/cameraman5x5.png)

Curvature regularization using 5x5 patches.

Installation
----------
MATLAB and a c++ compiler is needed. 
Setup matlab using mex -setup.
All files are compiled automatically.

The code has been tested on
* Matlab 2013a with GCC 4.8 on Ubuntu 13.10.
* Matlab 2013a with Visual Studio 2013 on Windows 7.

Getting started
----------
Two example are included in /examples
* binary_deconvolution_figure_12.m
* curvature_segmentation_figure_10.m

Included solvers
----------
The partial enumartion reformulation is optimized with

* [Vladimir Kolmogorovs' TRW-S implementation](http://pub.ist.ac.at/~vnk/papers/TRW-S.html).

The binary deconvolution problem is also optimized with

* [Vladimir Kolmogorovs' Roof duality implementation](http://pub.ist.ac.at/~vnk/software.html).

References
----------
1. Carl Olsson, Johannes Ulén, Yuri Boykov and Vladimir Kolmogorov 
 [Partial Enumeration and Curvature Regularization]
(http://www2.maths.lth.se/vision/publications/publications/view_paper.php?paper_id=584)
International Conference on Computer Vision. 2013.

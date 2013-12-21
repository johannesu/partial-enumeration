This is a MatLab/C++ library of our ICCV 2013 paper [1].

Installation
----------
You will need MATLAB and a c++ compiler. Setup matlab using mex -setup.
All files are compiled automatically.

The code has been tested on
* Matlab 2013a with GCC 4.8 on Ubuntu 13.10.
* Matlab 2013a with Visual Studio 2013 on Windows 7.

Usage
----------
Two example are included in /examples
* binary_deconvolution_figure_12.m
* curvature_segmentation_figure_10.m

Included solver
----------
The partial enumartion reformulation can be optimized with
many different optimizer, included is:

[Vladimir Kolmogorovs' TRW-S implementation](http://pub.ist.ac.at/~vnk/papers/TRW-S.html)

Roof duality is also included as a comparison in the binary deconvolution example.
[Vladimir Kolmogorovs' Roof duality implementation](http://pub.ist.ac.at/~vnk/software.html)

References
----------
1. [Partial Enumeration and Curvature Regularization](http://www2.maths.lth.se/vision/publications/publications/view_paper.php?paper_id=584)
Carl Olsson, Johannes Ulén, Yuri Boykov and Vladimir Kolmogorov
International Conference on Computer Vision. 2013.

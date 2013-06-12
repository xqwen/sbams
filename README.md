Companion Software for Bayesian Model Selection in Complex Linear Systems 
==========================================

This directory contains SBAMS, a package implementing Bayesian model selection in complex linear systems.

SBAMS is free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


Download
=============================================
The complete package can be downloaded from a single tar file in the download directory.


Source Code
=============================================
We provide R and C++ source code to compute approximate Bayes factors for Multivariate Linear Regression Models (MVLR), the R code can be found in src/R/mvlr.R, and the c++ code is in files src/mcmc/MVLR.cc and src/mcmc/MVLR.h

The directory src/mcmc/ also include c++ implemented Markov Chain Monte Carlo (MCMC) algorithm to perform Bayesian model selection for the MVLR models.  



Compilation and Installation
=============================================

The compilation from the c++ source code requires the GNU c++ compiler (g++), GNU make and GNU scientific library (GSL) pre-installed in the compiling machine. 

To compile the executable, run the following commands from the current directory:

     cd src/mcmc/
     make

Upon successful compilation, a binary executable "sbams_mvlr" should be produced.   


Documentation and Example Data
=============================================

A detailed documentation "sbams_mvlr.pdf" can be found in doc/ directory, we also include a simulated sample data set in the data/ directory.  


Citation
=============================================

Please cite the following for the usage of this software package:

Wen, X. "Bayesian Model Selection in Complex Linear Systems, as Illustrated in Genetic Association Studies", submit to Biometrics.

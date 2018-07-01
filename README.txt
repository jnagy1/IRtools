  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT AND OF THE ACCOMPANYING FILE LICENSE.txt.               *
  ***************************************************************************

   AUTHORS:

       Silvia Gazzola
       University of Bath, United Kingdom
       Email: s.gazzola@bath.ac.uk

       Per Christian Hansen
       Technical University of Denmark, Denmark
       Email: pcha@dtu.dk

       James G. Nagy
       Emory University, USA
       Email: jnagy@emory.edu

   REFERENCE:

       IR Tools: A MATLAB Package of Iterative Regularization
       Methods and Large-Scale Test Problems
       NUMERICAL ALGORITHMS, XX (XXXX), pp. XXX-XXX
       DOI: XXXXXXXXXXXXXXXXXXXXXXXXXXX.

   SOFTWARE REVISION DATE:

       V1.0, April 2018

   SOFTWARE LANGUAGE:

       MATLAB 9.3 (R2017b)


=====================================================================
SOFTWARE
=====================================================================

The IR Tools package for MATLAB provides iterative regularization methods
and test problems for large-scale linear inverse problems.

The software has been developed and tested using MATLAB version 9.3.
No other MathWorks products or toolboxes are required.

To obtain full functionality of this package, it is recommended to also
install the AIR Tools II package from:
   https://github.com/jakobsj/AIRToolsII


=====================================================================
HOW TO INSTALL AND CHECK THE INSTALLATION
=====================================================================

Please follow these steps:

- Extract the archive in a directory of your choice. This creates a
  directory called IRTools with all files of the toolbox.

- Start MATLAB.

- Add the directory IRTools containing the toolbox as well as all
  subdirectories to the MATLAB search path.
  - The easiest way to do this is by using the function
    IRtools_setup which updates the path permanently.
  - The alternative is to use the "Set Path" tool in the "Home" tab
    of the main MATLAB window, by clicking "Add with Subfolders..."
    and choosing the IRTools directory.

- The software installation can be checked by running example scripts
  with the generic name EX___.


=====================================================================
SOFTWARE UPDATES AND BUG FIXES
=====================================================================

In addition to the refereed version of the software published along
with the journal paper, the IRtools software will be maintained
in the GitHub code repository

  https://github.com/jnagy1/IRtools

Please check this location for software updates and bug fixes.


=====================================================================
PACKAGE
=====================================================================

The IR Tools package is organized into a main directory and subdirectories:

IRtools                    - contains this README.txt file, a LICENSE.txt
                             file, Contents.m which provides a detailed 
                             overview of all files in the package (can be 
                             listed from within MATLAB using "help IRTools"), 
                             and an installation function called IRtools_setup.m
IRtools/IRcodes            - contains all iterative solvers 
IRtools/PRcodes            - contains all test problems
IRtools/OPcodes            - contains operators needed for some of the test 
                             problems in PRcodes
IRtools/EXcodes            - contains some sample example scripts to 
                             illustrate how to setup and solve test problems
IRtools/Extra              - contains a number of auxiliary functions needed
                             in some of the main functions.
IRtools/Extra/PRblur_tools - contains some auxiliary functions needed
                             to handle debluring problems generated
                             with PRblur.
IRtools/Extra/test_data    - contains some data files for the deblurring
                             test problems.

See the file Contents.m for a list of the files contained in the package.


=====================================================================
EXAMPLE SCRIPTS
=====================================================================

We include a number of example scripts that can be executed for
illustration of the basic use and functionality of the package.
These scripts have the generic name EX___, and can be found in
IRtools/EXcodes

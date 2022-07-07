# Development
This document highlights current limitations/room for improvement to both the methodology as well as the software development. 
This is a good place to make note of things that need to be addressed by the creators and/or suggestions by the users. If you are a collaborator, please feel free to tackle any of the challenges listed here and submit a pull request.


## Methodology
These are potential areas of improvement or next steps in the development of the underlying methodology.

#### Sampling resonance parameters to create theoretical cross section
  - Spin group average parameter value sampling vs pulling from endf
    - currently pulling from endf, but need to correct PT distribution for very skinny resonances not being reported in endf
  - Feeding to sammy

#### Creating experimental data points from theoretical cross section
  - Appropriate energy spacing for experimental data points
  - Experimental noise sampling

## Software
These are potential areas of improvement or next steps in the development of this project as a software package. 

#### Modules
The syndat.sammy_interface.samtools_fmtpar() module runs into a strange issue maintaining the 11 character width

#### Packgaing and User installation
  - Look into putting dependencis into setup.cfg file rather than the current method of (pip installing all dependencies in requirements. file). This may be more robust.

#### Sphinx
  - Inlcude sphinx/latex extension that allows for mathematical expressions in documentation
  - Get example to render plot in html
  - In the examples for stochastic functions, use an independent RNG to set the seed rather than the depreciated global np.random.seed()
  



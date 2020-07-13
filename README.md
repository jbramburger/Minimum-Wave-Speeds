This repository contains the MATLAB files to reproduce the data and figures from "Minimum wave speeds in monostable reaction-diffusion equations: sharp bounds by polynomial optimization" by Jason J. Bramburger and David Goluskin (2020).

Computations use YALMIP version R20190425 to translate the SOS constraints into semidefinite programs and these are solved using Mosek version 9.0. These programs are publicly available at:

YALMIP: https://yalmip.github.io/download/
Mosek: https://www.mosek.com/downloads/

To reproduce data the scripts assume that YALMIP files are stored in the folder 'YALMIP-master' and Mosek files in the folder 'mosek'. 

In all scripts the parameter 'd' represents the degree of the auxiliary polynomial for either the surface or the volume method, while all other parameters specific to the differential equation are clearly delineated and have the same name as in the manuscript. Specifically, the scripts do the following:

- mFisher_Surface: Surface method applied to the modified Fisher equation, as given in Section 3.3. Implementation details are found in Section 3.1.

- mFisher_Volume_UpperBnd: Volume method for producing upper bounds on the minimum wave speed to the modified Fisher equation, as given in Section 3.3. Implementation details are found in Section 3.1.

- mFisher_Volume_LowerBnd: Volume method for producing lower bounds on the minimum wave speed to the modified Fisher equation, as given in Section 3.3. Implementation details are found in Section 3.2.

- chemotaxis_Surface: Surface method for producing upper bounds on the minimum wave speed to the chemotaxis model featured in Section 3.4. Implementation details are found in Section 3.1.

- chemotaxis_Volume: Volume method for producing lower bounds on the minimum wave speed to the chemotaxis model featured in Section 3.4. Implementation details are found in Section 3.2.

- autocatalytic_Volume_UpperBnd: Volume method for producing upper bounds on the minimum wave speed to the two-component autocatalysis reaction-diffusion equation featured in Section 4. Implementation details are found in Section 4.1.

- autocatalytic_Volume_LowerBnd: Volume method for producing lower bounds on the minimum wave speed to the two-component autocatalysis reaction-diffusion equation featured in Section 4. Implementation details are found in Section 4.1.



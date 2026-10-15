# Partition_ANW_UVVIS
---

The ANW partitioning model separates the non-water absorption coefficient, a<sub>nw</sub>, into the separate contributions of the phytoplankton absorption coefficient, a<sub>ph</sub>, and the non-phytoplankton absorption coefficient, a<sub>dg</sub>. The spectrum of a<sub>dg</sub> represents the sum of the non-algal (depigmented) particulate absorption coefficient, a<sub>d</sub>, and the chromophoric dissolved organic matter (CDOM) absorption coefficient, a<sub>g</sub>. The model partitions an input spectrum of a<sub>nw</sub> using a predefined library of spectral shape functions for a<sub>dg</sub> and initial inequality constraints associated with phytoplankon absorption band ratios to generate a pool of speculative solutions, which are then subjected to additional inequality or other constraints in order to define a pool of feasible solutions. The final step determines the model optimal solutions based on a statistical approach. This version of the model is configured for use in the UV-VIS spectral range of 350 to 700 nm and utilizes a spectral shape library specifically defined for this spectral range. 

A full description of the model can be found in the publication:
Stramski, D., and R. A. Reynolds. 2026. Model for partitioning the light absorption coefficient of seawater within the near-UV and visible spectral region into phytoplankton and non-phytoplankton components. Applied Optics, 65(30), 10467-10482. https://doi.org/10.1364/AO.610292

The current version of the model is v1.1. The model source code is in MATLAB file format. This README document provides information about the files within the Partition_ANW_UVVIS repository.

---

## Partition_ANW_UVVIS.m:
MATLAB function implementing the ANW partitioning model to partition input a<sub>nw</sub> into a<sub>dg</sub> and a<sub>ph</sub> using a predefined library of spectral shape functions of a<sub>dg</sub> in the UV-VIS along with a set of predefined inequality and other constraints. Returns output wavelength and solutions for a<sub>dg</sub> and a<sub>ph</sub>. This code includes multiple explanatory comments describing operation of the model.

## ADG_shapelibrary_UVVIS_498.mat:
MATLAB file containing the library of spectral shape functions defined in the UV-VIS spectral region (from 350 to 700 nm with a 1 nm spectral interval). See the supporting publication (in preparation) for further details about the formation of the library.

## TEST_Partition_ANW_UVVIS.m:
MATLAB script that executes Partition_ANW_UVVIS.m on a sample input spectrum of a<sub>nw</sub> defined from 350 to 700 nm in 1 nm intervals. The code utilizes the predefined UV-VIS spectral shape function library of a<sub>dg</sub> and the default values of eight inequality constraints to partition the input a<sub>nw</sub> into a<sub>dg</sub> and a<sub>ph</sub>.

## Partition_ANW_UVVIS_Test_Run.xls:
Microsoft Excel spreadsheet containing the input spectrum of a<sub>nw</sub> and resulting output data of a<sub>dg</sub> and a<sub>ph</sub> obtained from running the script TEST_Partition_ANW_UVVIS.m.

---

Contributors: Rick A. Reynolds, Dariusz Stramski, and Matthew Kehrli\
Contacts: Rick Reynolds (rreynolds@ucsd.edu), Dariusz Stramski (dstramski@ucsd.edu)
Ocean Optics Research Laboratory, Scripps Institution of Oceanography, University of California San Diego

---

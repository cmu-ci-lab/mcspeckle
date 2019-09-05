# SCMC: Speckle Covariance and rendering Monte-Carlo solver

A Monte Carlo rendering framework for the physically-accurate simulation of speckle patterns and statistics
arising from volumetric scattering of coherent waves.
This framework allows efficient computation of speckle mean and covariance as well as direct rendering of speckle fields,
for user defined light source directions and detector views. 
Only bulk parameters of the scattering volume are required as input, 
and there is no need to know the specific positions of wavelength-size scatterers in the medium.

Based on : ["A Monte Carlo Framework for Rendering Speckle Statistics in Scattering Media", Bar et.al. (2019)](https://arxiv.org/abs/1901.06931).

## Getting Started

These instructions will get your copy of the SCMC solver up and running on your local machine.

### Prerequisites

The SCMC solver requires MATLAB, available at: 

[https://www.mathworks.com/products/matlab.html](https://www.mathworks.com/products/matlab.html)

The solver was tested on MATLAB 2017b.

### Installing

Download the project code (available in github) to your local machine, and unzip all the files.

**For permanent installation**, add the SCMC folder and all its sub folders to the path, by

```
home -> environment -> set path -> add with Subfolders...
```

and select the SCMC folder. Then, restart MATLAB software.

**For temporal installation**, it is possible to add temporally the SCMS to the path, by adding the line

```
addpath(genpath(' '));
```

to your code, where the path to SCMC folder is inside the quotation marks.

## Running the code

### Interface

The scmc.m is the interface function of the solver.
Refer to the help of scmc.m and the examples folder for proper parameter specifications and usage.

### Examples

Several usage examples and applications are located in the examples subfolder. 
It is recommended to start with these basic examples.

### Assumptions
- The following scattering functions are supported: isotropic, Henyey-Greenstein (HG),
    and a user defined tabulated function.
- This code does not compute speckle mean. As the mean decays exponentially with the optical depth, in most cases it is  negligible. 

### Core code

The core algorithms are implemented in MCcov.m (evaluating speckle covariance) and MCfield.m (sampling a spackle field), which are located in code subfolder.

## Authors

* **Chen Bar** - *Department of Electrical Engeneering, Technion, Israel*
* **Marina Alterman** - *Department of Electrical Engeneering, Technion, Israel*
* [**Ioannis Gkioulekas**](http://www.cs.cmu.edu/~igkioule/) - *Robotics Institute, CMU, USA.*
* [**Anat Levin**](http://webee.technion.ac.il/people/anat.levin/) - *Department of Electrical Engeneering, Technion, Israel*

## License

This project is licensed under the ???

## Acknowledgments

* ?!?!?
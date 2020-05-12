![alt text](https://github.com/batflorez/MorphoSegger_v1/blob/master/Morphometrics_v2/Morphometrics_GUI/morphometrics_v2_icon.JPG?raw=true) 

# MorphoSegger v1.0

MorphoSegger is a Matlab package for imaging analysis of bacterial cells that combines two powerlful tools, [SuperSegger](https://github.com/wiggins-lab/SuperSegger/wiki) and [Morphometrics](https://simtk.org/projects/morphometrics). SuperSegger provides machine-learning based segmentation (with training GUI, foci and fluorescence analysus) and Morphometrics interpolates contours at subpixel resolution for robust cell size and morphology estimation. [Fiji](https://fiji.sc/) macros are integrated in Matlab through [MIJ](http://bigwww.epfl.ch/sage/soft/mij/) to expand the range of tools for pre-processing and analysis. This tool is most suited for analysis of bacteria in pads or CellAsic platforms. 

## Installation

#### Requirements:

  
  * Curve fitting Toolbox  
  * Statistics and Machine Learning Toolbox  
  * Global Optimization Toolbox  
  * Parallel Computing Toolbox  
  * Running installation of [Fiji](https://fiji.sc/) and [MIJ](http://bigwww.epfl.ch/sage/soft/mij/)  
  * Adapted code of SuperSegger and Morphometrics (see below for a list of modifcations).  
**Note:** set Matlab memory settings: *Preferences > General > Java Heap Memory* 


#### Setting up MIJ:
The easiest way to install MIJ is via the repository. Click *Update > Manage update sites* on Fiji and check the repository [BIG-EPFL](https://sites.imagej.net/BIG-EPFL/). It might be useful to add the [ImageJ-Matlab](https://sites.imagej.net/MATLAB/) repository as well. Finally, add the path for Fiji scripts in Matlab: `addpath '/Applications/Fiji.app/scripts'`. Type `Miji` in Matlab to make sure it loads properly.
  
#### Modifications to SuperSegger and Morphometrics:





Download the tool:

```
git clone https://github.com/batflorez/MorphoSegger.git
```

## Usage

Prepare a pipeline script for each experiment with the conditions and analysis options needed. Then run in Matlab command window:
```
procesExp
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate.

## License
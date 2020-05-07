# MorphoSegger v1.0

MorphoSegger is a Matlab package for imaging analysis of bacterial cells that combines two powerlful tools, [SuperSegger](https://github.com/wiggins-lab/SuperSegger/wiki) and [Morphometrics](https://simtk.org/projects/morphometrics). It uses ImageJ macros run in Matlab with [MIJ](http://bigwww.epfl.ch/sage/soft/mij/) for image pre-processing.

## Installation

### Requirements:
It requres the following Matlab toolboxes:

Image Processing Toolbox  
Deep Learning Toolbox  
Curve fitting Toolbox  
Statistics and Machine Learning Toolbox  
Global Optimization Toolbox  
Parallel Computing Toolbox  
Running installation of [Fiji](https://fiji.sc/) and [MIJ](http://bigwww.epfl.ch/sage/soft/mij/)

Dowonload the tool:

```bash
git clone https://github.com/batflorez/MorphoSegger_v1.git
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
[MIT](https://choosealicense.com/licenses/mit/)
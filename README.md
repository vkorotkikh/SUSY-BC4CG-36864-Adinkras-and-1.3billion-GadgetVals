# BC4 Adinkra Calculation Code
Description: This repository contains the python scripts thatenable the calculation
of all all unique 36,864 ordered BC4-based Adinkras with four colors, four
open-nodes and four closed nodes and then calculates the Vij holoraumy matrices
and the 1,358,954,496 Matrix elements of the Gadget.
My original, basic version of the code used for calculations in the paper is
found under original-code branch.

## DEV STATUS
```
Currently upgrading master branch w/ enhancements which may (T/F)? break current code.
```


## Important details
Developed and tested using Python 3 (using Anaconda Python 3.6 as of README update).
```
$ Requires numpy library
$ Developer/creator tested on Mac OSX. Untested for Windows usage.
```

## Getting Started
To execute this tool, run the following command in Terminal
```
$ python -u adinkra_tetrad_calc.py
```
or if you have both python 2 and 3, specify python 3 version.
```
$ python3 adinkra_tetrad_calc.py
```


## General Overview
### Primary script
adinkra_tetrad_calc.py - Calculates the 36,864 ordered BC4-based
adinkras with 4 colours, 4 open and 4 closed nodes. In this case each Adinkra is
a tetrad of 4 L matrices.
### Code for Gadget calculation
vij_holoraumy_calc.py - Calculates the Vij holoramy matrices from each Adinkra.
There's six Vij matrices per each Adinkra. The script proceeds to calculate
the billion plus matrix elements of the Gadget.

## Author

* **Vadim K.** - **

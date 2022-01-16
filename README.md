# Overview

This project contains the source for building force-field components for
the phenix crystallographic refinement program.  It contains"msander",
a "modern" version of parts of sander, plus other pieces of AmberTools
needed for basic building and simulations of biomolecules.  Tools inlcuded
are:
```
   addles  antechamber  tleap  msander  parmed  sqm
```
Also included are the API's to msander, and various X-ray-related utilities

# Building the code

*Conda build
```
   conda build [ --python x.x ] recipe 
      (note: you should have conda-forge at the top of your channel
      list in ~/.condarc.  You should also have done a "conda install
      conda-forge-pinning" in your conda build environment.
```

*Non-conda build  (MacOSX, Linux, probably WSL):
```
   ./configure --help   #  then choose the options you want
   make install
   make test
```

# License
This project is generally licensed under the GNU (Lesser) General Public 
License, version 3 (GPL/LGPL v3).  Some components use different, but 
compatible, open source licenses.  See the LICENSE file for more information.


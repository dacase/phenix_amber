# Overview

This directory tree contains "msander", a "modern" version of parts of
sander, plus other pieces of AmberTools needed for basic building and
simulations of biomolecules.  Tools inlcuded are:
```
   addles  antechamber  tleap  msander  parmed  sqm
```
Also included are the API's to msander, and various X-ray-related utilities

# Warning

This is a work in progress, and may not always be in a stable
state.  I may not be able to respond to requests for support.
The documentation is here:

    https://ambermd.org/doc12/Amber21.pdf

This code is probably only useful to those who are already familiar with
AmberTools, and there are some small differences that are not yet
documented.  You can look in the test directory for examples of input files,
or send email to dacase1@gmail.com if you want to participate in development.

# Design goals

* This project began as a fork of the `sander` code in `AmberTools`.  
It tries to (greatly) simplify the code base, choosing the best and 
most useful parts of the code, and to serve as a test bed for how 
modern Fortran coding techniques can be used.  Key application areas 
are expected to be in structure refinements using NMR, cryoEM or 
Xray diffraction information.  This version has a fair amount of OpenMP
support, especially for Xray and 3D-RISM calculations.

* This project is also the source for the phenix_amber project, which explains
why non-sander things like tleap, antechamber, etc. are here.

* Since this code is based on sander, tons of people have been involved in its
creation over the years.  See https://ambermd.org/contributors.html for more
information, although even that gets out of date.

# Key differences in functionality versus sander

* Some pieces are missing from the sander program in AmberTools:

  * Things that should be easy to re-introduce later: emil, sebomd, pbsa, APBS

  * Things are are problably gone for good, but which don't represent the best
current practice: Path-integral methods, thermostats that don't follow
the "middle" scheme, Berendsen barostat

  * Things that might be useful, but really complicate the code: evb
potentials, some parts of adaptive QM/MM, nudged elastic band, constant pH
and constant redox potential simulations.

  * Non-periodic 3D-RISM has been removed for now, in an attempt to get the
simplest possible RISM code, perhaps as a basis for future GPU work.

* Key pieces of code that are still there, and being emphasized:

  * Periodic and non-periodic simulations, with all of Amber's GB models

  * QM/MM, including hooks to external codes

  * NMR, cryoEM and Xray restraints (including quite a bit of new code)

  * Thermodynamic integration and non-equilibrium sampling methods

  * Replica exchange capabilities, except for constant pH and redox potential
simulations

# Building the code

*Conda build (serial mode only, no MPI):
```
   conda build [ --python x.x ] recipe 
      (note: you should have conda-forge at the top of your channel
      list in ~/.condarc.  This build creates the conda package used
      in the phenix_amber project, but the package could also have other uses.)
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


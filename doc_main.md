Main Page                                    {#mainpage}
=========

# Introduction #

Software `argon39` handles the calculation of nucleogenic production rates of neutrons, <sup>21</sup>Ne and <sup>39</sup>Ar in underground geological settings. 

`argon39` is not a standalone tool. It requires interaction with `TALYS` ([www.talys.eu](http://www.talys.eu)) and `MCNP6` ([mcnp.lanl.gov](http://mcnp.lanl.gov)) codes.

`argon39` was used for calculations presented in

*  O. Šrámek, L. Stevens, W. F. McDonough, S. Mukhopadhyay, and J. Peterson: “Subterranean production of neutrons, <sup>39</sup>Ar and <sup>21</sup>Ne: Rates and uncertainties.” In preparation.

Written in Fortran 90/95 with some Fortran 2003/2008 features.


# Quick Start #

Download or clone the repository on your machine. Go to the main directory (`argon39/`). Makefile is included, default Fortran compiler is set to `gfortran`. To compile the code, in a unix shell type

    $ make 

Four executables should be created: `argon39.x`, `mcnp6extract.x`, `mcnp6prep.x`, `talysprep.x`.

Change to directory `example/` which contains an input file `argon39_example.inp`. Run the code by

    $ ../argon39.x argon39_example.inp

A summary output file `argon39_example.out` will be created + many other output files.

Please see [Detailed Description](@ref detail) for more.


# Project Info #

@date     April, 2013 -- September, 2015

@author   [Ondřej Šrámek](http://www.ondrejsramek.net), Charles University in Prague, Department of Geophysics, ondrej.sramek@gmail.com

Collaboration with 

* [Bill McDonough](https://www.geol.umd.edu/~mcdonoug/) (Univ. Maryland)
* [Sujoy Mukhopadhyay](http://mygeologypage.ucdavis.edu/sujoy/) (UC Davis)
* [Jerry Peterson](http://phys.colorado.edu/people/peterson-jerry) (CU Boulder)
* [Lauren Stevens](http://www2.chem.umd.edu/groups/beichhorn/Eichhorn_Research_Group/People/Pages/Lauren.html) (Univ. Maryland)

Built and tested on Mac OS X 10.9 with [GNU Fortran](http://gcc.gnu.org/fortran/) (`gfortran`) gcc version 4.9.3 and 5.2.0 and [GNU make](http://www.gnu.org/software/make/) 3.81. Documentation generated from annotated F95 sources and Markdown text files using [Doxygen](http://doxygen.org) version 1.8.9.1. Version control provided by [Mercurial](http://mercurial.selenic.com/) version 3.5 and Git version 2.5.0.


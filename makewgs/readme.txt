

MAKE WING OR BODY FILES IN LAWGS FORMAT                      /makewgs/readme.txt

The files for this program are in the directory makewgs 
  readme.txt     general description
  input.txt      description of the input file for makewgs
  makewgs.f90    the complete source code in modern Fortran


Sample cases for this program are:
  delta1.mak     input data for simple delta wing
  sh.mak         input data for Sears-Haack body
  delta1.wgs     LaWGS file produced by delta1.mak
  sh.wgs         LaWGS file produced by sh.mak

The reference documents for this program may be accessed
from the web page https://www.pdas.com/makewgsrefs.html. 

To compile this program for your machine, use the command
   gfortran  makewgs.f90 -o makewgs.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

This program simply asks for the name of the input file. This must be 
a file containing commands described in input.txt.
Upon successfuly reading the input data, the program produces a file
called make.wgs that may be used as input to programs requiring a
wireframe geometry file.

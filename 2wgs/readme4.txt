

CONVERT SHABP INPUT FILE TO LAWGS FORMAT                  /2wgs/readme4.txt

The files for this program are in the directory  2wgs
  readme4.txt     General description
  hab2wgs.f90     the complete source code (still 'In Progress').

Sample cases for this program are:
  tnd6480.mk5     input  data for wing-body configuration of NASA TN D-6480
  tp2608.mk5      input  data for wing-body configuration of NASA TP 2608
  tnd6480.wgs     output data for wing-body configuration of NASA TN D-6480
  tp2608.wgs      output data for wing-body configuration of NASA TP 2608

The reference documents for the tnd6480 and tp2608 configurations may be accessed
from the web page https://www.pdas.com/2wgsrefs.html. The NASA reference for
LaWgs may also be found there.
 
This program asks for the name of the input file. This must be a file written 
in the input format for the Hypersonic Arbitrary Body Program.
After reading the input data, the program produces a file called hab.wgs that may
be used as input to 3view, hlp or wgs2wrl.

To compile this program for your computer, use the command
   gfortran hab2wgs.f90 -o hab2wgs.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

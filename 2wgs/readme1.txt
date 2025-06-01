

CONVERT WINGBODY INPUT FILE TO LAWGS FORMAT               /2wgs/readme1.txt

The files for this program are in the directory 2wgs 
  readme1.txt     general description
  wb2wgs.f90      the complete source code

Sample cases for this program are:
  tnd4211.inp     input  data for wing-body configuration of NASA TN D-4211
  tnd7505.inp     input  data for wing-body configuration of NASA TN D-7505
  tnd4211.wgs     output data for wing-body configuration of NASA TN D-4211
  tnd7505.wgs     output data for wing-body configuration of NASA TN D-7505

 
This program asks for the name of the input file. This must be 
a file written in the input format to WingBody. After reading 
the input data, the program produces a file called wb.wgs that may be 
used as input to 3view, hlp or wgs2wrl.

To compile this program for your machine, use the command
   gfortran  wb2wgs.f90 -o wb2wgs.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

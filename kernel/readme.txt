

A STEADY AND OSCILLATORY KERNEL FUNCTION METHOD            /kernel/readme.txt
FOR INTERFERING SURFACES IN SUBSONIC, TRANSONIC AND SUPERSONIC FLOW


The files for this program are in the directory kernel 
  readme.txt      general description
  input.txt       input description (Appendix H from NASA CR-144895)
  kernel.f90      the complete source code in modern Fortran
  lar12524.txt    the original program description from COSMIC
  original.src    the original copy of the source code (from COSMIC)

Sample cases for this program are:
  case1.inp       input data
  case1.out       output data for case1.inp

The reference documents for this program may be accessed
from the web page https://www.pdas.com/kernelrefs.html. 

To compile this program for your machine, use the command
   gfortran  kernel.f90 -o kernel.exe
Linux and Macintosh users may prefer to omit the .exe on the output file name.

This program first asks for the name of the input file. This must
be a file written to conform to the input instructions.  After calculating 
the solution, the program produces a file called kernel.out that contains a 
wealth of information concerning the flow problem to be solved.

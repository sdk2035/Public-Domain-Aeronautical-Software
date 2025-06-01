

ABAXI - TRANSIENT RESPONSE OF ABLATING AXISYMMETRIC            /abaxi/readme.txt
 BODIES INCLUDING THE EFFECTS OF SHAPE CHANGE

The files for this program are in the directory abaxi 
  readme.txt      general description
  original.src    the original source code (from COSMIC)
  abaxi.f90       the complete source code in modern Fortran
  lar11049.txt    the original program description (from COSMIC)

Sample cases for this program are:
  tmcase1.inp     input data for test case 1 (from TM X-2375)
  tmcase1.out     output data for tmcase1.inp
  tmcase2.inp     input data for test case 2
  tmcase2.out     output data for tmcase2.inp
  tmcase3.inp     input data for test case 3
  tmcase3.out     output data for tmcase3.inp

The NASA reference documents TM X-2375 and TN D-6220 may be accessed
from the web page https://www.pdas.com/abaxirefs.html and the description
of the input to ABAXI will be found in TM X-2375.

This program first asks for the name of the input file. This must
be a file written to conform to the input definition.  After calculating 
the solution, the program produces a file called abaxi.out that contains a 
wealth of information concerning the flow problem to be solved.
The lines in abaxi.out are 132 characters in length, so you need to set
your printer to landscape mode, 8 point characters to keep the output
on the page.

To compile this program for your machine, use the command
   gfortran  abaxi.f90 -o abaxi.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.



MASSPROP- MASS PROPERTIES OF A RIGID STRUCTURE              /massprop/readme.txt

The files for this program are in the directory massprop
  readme.txt      general description
  input.txt       guide for writing the input file in PDAS format
  massprop.f90    the complete source code in modern Fortran
  lar12454.txt    the original program description from COSMIC
  original.src    the original copy of the source code (from COSMIC)

Sample cases for this program are:
  example1.inp    input data for example 1
  example1.out    output data for example1.inp
  example1.dbg    debug data for example1.inp 
  example2.inp    input data for example 2
  example2.out    output data for example2.inp
  example3.inp    input data for example 3
  example3.out    output data for example3.inp

The reference documents for this program may be accessed
from the web page https://www.pdas.com/massproprefs.html. 

You may execute the program by typing the command massprop followed
by the filename of the input file. If you omit the filename
(or misspell it), the program prompts you for the name of the input file. 
This must be a file written to conform to the PDAS format as 
described in the file input.txt.  This has been modified from the original
format described in the reference document. After the calculations are
complete, the program produces a file called massprop.out that contains a 
summary of the results. The program also produces a file named massprop.dbg
that displays some intermediate results for each component that can be
useful in making a detailed examination of the case.


To compile this program for your machine, use the command
   gfortran  massprop.f90 -o massprop.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

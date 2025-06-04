
SONIC BOOM                                                      /boom/readme.txt
   
The files for this program are in the directory boom 
  readme.txt      general description
  original.src    the original source as scanned from NASA CR-157
  boom.f90        the complete source code in modern Fortran
  case1.inp       input for sample case 1
  case1.out       the output from case1.inp
  case2.inp       input for sample case 2
  case2.out       the output from case2.inp
                
The reference documents for this program may be accessed from the web page 
https://www.pdas.com/boomrefs.html and the description of the input file may 
be accessed from NASA CR-157.
 
This program asks for the name of the input file. This must be 
a file written in the input format to Boom. After reading 
the input data, the program produces a file called boom.out containing
the program output for printing.

To recompile this program for your computer, use the command
   gfortran  boom.f90 -o boom.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

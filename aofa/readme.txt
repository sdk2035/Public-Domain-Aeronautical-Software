

AOFA- THREE-DIMENSIONAL SUPERSONIC VISCOUS FLOW                 /aofa/readme.txt

The files for this program are in the directory aofa 

  readme.txt      general description
  original.src    the original copy of the source code (from COSMIC)
  arc11087.txt    the original program description from COSMIC
  aofa.f90        the source code converted to modern Fortran


This program is a "work in progress" and is not ready for general release.
I have included it so those who have a special interest may see
the original code plus my modifications to date.

To compile this program for your computer, use the command
   gfortran aofa.f90 -o aofa.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

The program asks for the name of the input file. This must be 
a file written in the input format to aofa. After reading 
the input data, the program produces a file called aofa.out containing
the program output for printing.


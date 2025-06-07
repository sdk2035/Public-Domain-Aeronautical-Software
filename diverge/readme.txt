

AEROELASTIC DIVERGENCE CHARACTERISTICS OF UNGUIDED,         /diverge/readme.txt
   SLENDER BODY, MULTI-STAGE LAUNCH VEHICLES

The files for this program are in the directory diverge 
  readme.txt      this file of general description
  original.src    the original copy of the source code (from COSMIC)
  diverge.f90     the complete source code in modern Fortran

One sample case for this program has been supplied:
  case1.inp       input data
  case1.out       output data for case1.inp
  case2.inp       same data, slightly rearranged
  case2.out       output data for case2.inp

 
This program asks for the name of the input file. This must be 
a file written in the input format to Diverge. Unfortunately, there is
no reference to define the data and format. The program produces a file 
called diverge.out that may be printed. This program may be very useful,
but I can't tell what it does.

To compile this program for your machine, use the command
   gfortran  diverge.f90 -o diverge.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

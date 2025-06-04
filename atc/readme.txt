

ANALYSIS OF AIRCRAFT MOTIONS                                     /atc/readme.txt

The files for this program are in the directory atc
  readme.txt      general description
  input.txt       guide to preparing input for atc
  atc.f90         the complete source code in modern Fortran
  original.src    the original copy of the source code (from COSMIC)
  arc11132.txt    the original program description from COSMIC    

Files for this sample case are:
  case1.inp       general input data
  radar.inp       sample radar data
  wind.inp        wind and temperature at altitude
  plot.inp        plotting instructions
  case1.out       output data for case1.inp and associated files

The program will ask the user for the names of the input data files. 
After reading the input data, the program produces atc.out with the 
tabulated data and printer plots.

To compile this program for your machine, use the command
   gfortran  atc.f90 -o atc.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

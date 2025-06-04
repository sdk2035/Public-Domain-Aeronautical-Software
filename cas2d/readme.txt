

CAS2D- NONROTATING BLADE-TO-BLADE, STEADY,                    /cas2d/readme.txt
 POTENTIAL TRANSONIC CASCADE FLOW ANALYSIS CODE


The files for this program are in the directory cas2d 


  readme.txt      this file of general description
  original.src    the original copy of the source code (from COSMIC)
  lew13854.txt    the original program description from COSMIC
  cas2d.f         F77 version of the original source
  case1.inp       a test case 
  case1.out       output from running case1.inp
The reference documents for this program may be referenced from the web
page at http://www.pdas.com/cas2drefs.html  

To compile this program for your computer, use the command
   gfortran cas2d.f90 -o cas2d.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

To execute this program using infile as input and producing outfile
as output, use the command

   cas2d < infile > outfile 

This program is a "work in progress" and is not ready for general release.
I have included it so those who have a special interest may see
the original code plus my modifications to date.
 

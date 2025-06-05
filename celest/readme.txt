

CELEST - TRANSFORMATION OF CELESTIAL COORDINATES              /celest/readme.txt

The files for this program are in the directory celest

  readme.txt      this file of general description
  celest.f90      the source code for the module of transformation equations
  tables.f90      a sample program that produces tables that are a subset
                     of the extensive ones in the document
  table1.txt      a table of equatorial to ecliptic transformation
  table2.txt      a table of ecliptic to equatorial transformation
  table3.txt      a table of ecliptic to galactic transformation
  table4.txt      a table of equatorial to galactic transformation
  table5.txt      a table of galactic to equatorial transformation

 
If you compile and run the program tables, it will produce the five text
files of transformation tables shown above. The file tables.f90 will
include the computational module in celest.f90.

  Windows:       gfortran tables.f90 -o tables.exe
  Mac or Linux:  gfortran tables.f90 -o tables

The CELEST procedures relate three basic frames of reference for finding position.
1. Equatorial (geocentric) system. The equatorial plane, determined
by the plane through the earth's equator, describes a great circle on
the celestial sphere. 
2. Ecliptic (heliocentric) system. The ecliptic plane is determined by the
path of earth's orbit around the sun in the solar system. The great circles
of the equatorial and ecliptic systems intersect at the First Point of Aries. 
This is considered to be the (0,0) point for these two systems. 
3. Galactic system. The galactic plane is formed by the Milky Way. 
The Milky Way is a relatively flat galaxy, so for our purposes, 
it can be adequately represented by a plane. 

Some Fortran subroutines, required to give the coordinate transformations
among the equatorial, ecliptic, and galactic systems, have been defined and coded.
A brief review of spherical trigonometry is included in the reference document
NASA Technical Memorandum 53943. (See https://www.pdas.com/celestrefs.html for a
link to download this document).
These coordinate transformations were developed and programmed
for use in the determination of the ecliptic (Zodiacal Light) with respect to the
other systems for analysis of the S-073/T-027 AAP experiment data analysis.

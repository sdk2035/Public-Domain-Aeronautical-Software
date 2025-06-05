

CONPLOT - A GENERAL ALGORITHM FOR THE                        /conplot/readme.txt
   CONSTRUCTION OF CONTOUR PLOTS               

The files for this program are in the directory conplot 
  readme.txt      general information
  arc11441.txt    original COSMIC description file
  conplot.f90     source code
  original.src    the original source code from COSMIC
  test1.f90       a sample program using a function coded in Fortran
  test2.f90       a sample program that reads data
  test2dat.txt    some example data for test2
  test2con.txt    some example contours for the data in test2dat.txt

To use this program, make a directory on your hard disk and copy these 
files to that directory.

The principal items of code are the modules Smooth and Contour.
Conplot is not a program, but a module that you use in your programs.
I have supplied two programs that let you practice.

test1 plots contours using data points that you compute in the
main program. There is no input data. It makes a file called
test.gnu that may be used with gnuplot.

test2 plots contours using data points that are read from a file.
Each line of the input data contains three numbers - x, y, f(x,y).
The program reads until there is no more data and then inquires
about a file of contour values. This file contains one value per record
containing the contour level for each contour line to be plotted.
There is a sample case called test2dat.txt that represents the pressure
coefficients on a delta wing in supersonic flow. It is taken from
NASA TN D-1264 by Harry Carlson. The file called test2con.txt has some
contour values that may be used. You can experiment with smoothing
and see the change in contours. To begin, try saying that you want 
smoothing and then try iexp=3 and jexp=4. I am not really sure what 
these mean, but you might figure it out from the source code.

Test2 produces a file called test.gnu that may be used with gnuplot.
For example,

   gnuplot>plot 'test.gnu' with lines

To avoid distortion of the image use the plotting command
   gnuplot> set size ratio -1


DESCRIPTION (from COSMIC) -
The graphical presentation of experimentally or theoretically generated data
sets frequently involves the construction of contour plots. A general
computer algorithm has been developed for the construction of contour plots.
The algorithm provides for efficient and accurate contouring with a modular
approach, which allows flexibility in modifying the algorithm for special
applications. The algorithm accepts as input data values at a set of points
irregularly distributed over a plane. The algorithm is based on an
interpolation scheme in which the points in the plane are connected by
straight line segments to form a set of triangles. In general, the data is
smoothed using a least-squares-error fit of the data to a bivariate
polynomial. To construct the contours, interpolation along the edges of the
triangles is performed, using the bivariable polynomial if data smoothing was
performed. Once the contour points have been located, the contour may be
drawn.

COMMENTS -
There are many other algorithms available for computation of contour lines.
There are several in the Transactions for Mathematical Software and 
elsewhere. This routine was deemed worthy of inclusion in the NASA COSMIC
collection, but I cannot swear to its status among contour generators.

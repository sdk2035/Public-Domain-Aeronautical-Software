

LINE INTERPOLATER - AN INTERPOLATION                          /linIntrp/readme.txt
                          AID FOR AIRPLANE LAYOUT       

The files for this program are in the directory linIntrp 
  readme.txt      general information
  int.dfm         the main form
  int.pas         the code for unit main 
  interp.dpr      the Delphi project file
  interp.exe      the executable file
 
In order to recompile this program, you will need to have Delphi installed
on your computer. You can use the executable file interp.exe without Delphi.

Line Interpolater is a simple tool that helps you compute all the points
you need for input of an airplane geometrical description to an engineering
analysis program. One frequently has a line in 3D defined by two points.
But additional points on this line are needed for input. You may know one
of the coordinates, but need the other two. Sometimes, you know the fractional
distance of the point between the two defining points.

This tool helps you get the numbers you need without laborious calculation.
On the left side of the screen are six boxes for entering the (x,y,z) 
coordinates of two distinct points that serve to define the line. On the 
right side of the screen are four input boxes that are set up as a radio
group - that is, only one is actually active. When you click the Compute
button (or press Enter), the point on the line corresponding to the active
coordinate is used to determine the remaining unknown coordinates.
 
I wrote this program as a new project and is not based on a program
from a national laboratory. The whole program, source code and all is
declared to be open source - public domain.
